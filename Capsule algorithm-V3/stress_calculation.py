'''
stress_calculation.py

包含计算由注浆引起的附加水平应力的函数。
主要基于论文公式 (7) (Mindlin 解的镜像法部分)，并应用修正简化法的修正系数。
'''
import numpy as np
import constants as const # type: ignore
from constants import GroutPoint # Import GroutPoint
from typing import List # Import List

def _calculate_sigma_x_from_point_source(x_obs, y_obs, z_obs, x_src, y_src, z_src, E, nu, alpha_strength, omit_complex_surface_terms: bool = False):
    '''
    根据论文公式 (7) 计算单个点状扩张源在观察点产生的水平附加应力 σ'_x(1-2)。

    参数:
        x_obs, y_obs, z_obs: 观察点坐标 (m)
        x_src, y_src, z_src: 源点坐标 (m) (z_src 是深度, >0)
        E: 土壤弹性模量 (Pa)
        nu: 土壤泊松比
        alpha_strength: 源的扩张强度 a^3(i-i0) (m^3)
        omit_complex_surface_terms: 是否省略公式(7)中与z相关的镜像项(term3_z_r2, term4_z0_r2)

    返回:
        sigma_x_contrib: 该点源贡献的 σ'_x(1-2) (Pa)
    '''
    if alpha_strength == 0:
        return 0.0

    # 相对坐标
    dx = x_obs - x_src
    dy = y_obs - y_src
    # dz_real = z_obs - z_src # z_obs 是观察点深度, z_src 是源深度

    # 到实际源的距离 r1
    # r1^2 = dx^2 + dy^2 + (z_obs - z_src)^2
    # r1 = np.sqrt(dx**2 + dy**2 + (z_obs - z_src)**2)
    # Ensure r1 is not zero to avoid division by zero
    # Smallest distance can be if obs point is at source, but source shouldn't be at tunnel itself.
    # Let's use a small epsilon for stability if r1 or r2 are extremely small, though physically they shouldn't be zero.
    epsilon = 1e-9

    # Terminology from paper for r1, r2 for clarity with Eq (7)
    # x, y, z in Eq (7) are field point (obs)
    # x0, y0, z0 in Eq (7) are source point (src)
    # So, (x-x0) in paper is dx, (y-y0) is dy, (z-z0) is (z_obs-z_src)

    r1_sq = dx**2 + dy**2 + (z_obs - z_src)**2
    if r1_sq < epsilon**2:
        r1 = epsilon
    else:
        r1 = np.sqrt(r1_sq)

    # 到镜像源的距离 r2 (镜像源在 -z_src)
    # r2^2 = dx^2 + dy^2 + (z_obs + z_src)^2
    # r2 = np.sqrt(dx**2 + dy**2 + (z_obs + z_src)**2)
    r2_sq = dx**2 + dy**2 + (z_obs + z_src)**2
    if r2_sq < epsilon**2: # Should always be >0 if z_obs > 0 and z_src > 0
        r2 = epsilon
    else:
        r2 = np.sqrt(r2_sq)

    # 防止 r1 或 r2 过小导致数值问题
    if r1 < epsilon or r2 < epsilon:
        # This case should ideally not be hit if source is reasonably far from tunnel
        return 0.0

    common_factor = (alpha_strength * E) / ((1 + nu) * (1 - 2 * nu) * 2 * np.pi)
    if common_factor == 0: # can happen if E=0 or alpha_strength=0
        return 0.0
        
    term1_r1 = 0.0
    term2_r2 = 0.0
    term3_z_r2 = 0.0
    term4_z0_r2 = 0.0

    # --- Contribution from real source (terms with r1) ---
    # {(1-2μ)((x-x0)^2/r1^5 - 1/r1^3) + 3(x-x0)^2/r1^5}
    # (x-x0) is dx
    r1_cubed = r1**3
    r1_fifth = r1**5
    if r1_fifth != 0: # should be covered by r1 < epsilon
        term1_r1 = (1 - 2 * nu) * (dx**2 / r1_fifth - 1 / r1_cubed) + 3 * dx**2 / r1_fifth
    
    # --- Contribution from image source (terms with r2) ---
    if not omit_complex_surface_terms:
        # Part 1: {(1-2μ)((x-x0)^2/r2^5 - 1/r2^3) + 3(x-x0)^2/r2^5}
        r2_cubed = r2**3
        r2_fifth = r2**5
        if r2_fifth != 0: # should be covered by r2 < epsilon
            term2_r2 = (1 - 2 * nu) * (dx**2 / r2_fifth - 1 / r2_cubed) + 3 * dx**2 / r2_fifth

            # Part 2: {(z/r2^3)((1-2μ)-3(x-x0)^2/r2^2)}
            # Here, z is z_obs (field point depth)
            term3_z_r2 = (z_obs / r2_cubed) * ( (1 - 2 * nu) - 3 * dx**2 / r2_sq )

            # Part 3: - {(z0/r2^3)((1-2μ)-3(x-x0)^2/r2^2 - (6*z*z0/r2^2)*(1 - (x-x0)^2/r2^2 - (z+z0)^2/r2^2))}
            # Here, z0 is z_src (source depth), z is z_obs (field point depth)
            # The term (1 - (x-x0)^2/r2^2 - (z+z0)^2/r2^2) simplifies to (y-y0)^2/r2_sq
            val_in_last_bracket = (1 - 2 * nu) - (3 * dx**2 / r2_sq) - \
                                  (6 * z_obs * z_src / r2_sq) * (dy**2 / r2_sq)
            term4_z0_r2 = - (z_src / r2_cubed) * val_in_last_bracket
    else: # omit_complex_surface_terms is True
        # Only include the simpler r2 term (term2_r2) if not omitting all surface terms.
        # For this specific simplification, we are omitting term3 and term4.
        # Let's re-evaluate if term2_r2 should also be omitted or if the request was only for terms 3 and 4.
        # Based on "omitting the more complex surface correction terms related to z and z0",
        # term2_r2 (Term B) is the basic image source term and should likely be kept.
        # The prompt was "省略了图中这一项", if "这一项" is specific, it might mean only one of them.
        # My interpretation was "term3_z_r2 and term4_z0_r2". Let's stick to that.
        # So if omit_complex_surface_terms is True, term3 and term4 are 0. term1 and term2 are still calculated.
        # Re-thinking: The logic of "omitting complex surface terms" should only zero out term3 and term4.
        # term2 (the basic image term) should still be calculated.
        r2_cubed = r2**3
        r2_fifth = r2**5
        if r2_fifth != 0:
            term2_r2 = (1 - 2 * nu) * (dx**2 / r2_fifth - 1 / r2_cubed) + 3 * dx**2 / r2_fifth
        # term3_z_r2 remains 0.0
        # term4_z0_r2 remains 0.0
        
    sigma_x_contrib = common_factor * (term1_r1 + term2_r2 + term3_z_r2 + term4_z0_r2)
    return sigma_x_contrib

def calculate_total_sigma_x_profile(
    grouting_points_profile: List[GroutPoint], # New parameter for list of grout points
    y_coords_tunnel,
    tunnel_x_center = const.X_TUNNEL_CENTER,
    tunnel_z_center = const.Z_TUNNEL_CENTER,
    calculation_method: str = "revised_simplified",  # 新增：计算方式选择
    override_eta_factor: float | None = None
):
    '''
    计算沿隧道轴线的由多个注浆点引起的总附加水平应力剖面 σ'_x。
    每个注浆点可以独立设置为均匀或非均匀线源。

    参数:
        grouting_points_profile: 注浆点参数列表 (List[GroutPoint])
        y_coords_tunnel: 隧道沿线计算应力的y坐标数组 (m)
        tunnel_x_center: 隧道中心x坐标 (m)
        tunnel_z_center: 隧道中心z坐标 (深度, m)
        calculation_method: 计算方式选择
            - "original": 原始计算（包括公式8的完整Mindlin解，η=1.0）
            - "simplified": 简化计算（省略复杂表面项，η=1.0）
            - "revised_simplified": 修正简化计算（省略复杂表面项，η=1.25）
        override_eta_factor: 如果提供，则覆盖自动选择的eta因子

    返回:
        sigma_x_profile: 沿y_coords_tunnel的附加水平应力数组 (Pa)
    '''
    
    # 根据计算方式确定参数
    if calculation_method == "original":
        omit_complex_surface_terms = False  # 完整Mindlin解
        eta_factor = 1.0
    elif calculation_method == "simplified":
        omit_complex_surface_terms = True   # 省略复杂表面项
        eta_factor = 1.0
    elif calculation_method == "revised_simplified":
        omit_complex_surface_terms = True   # 省略复杂表面项
        eta_factor = const.ETA_CORRECTION  # 使用constants中的修正系数
    else:
        raise ValueError(f"Unknown calculation_method: {calculation_method}. Must be 'original', 'simplified', or 'revised_simplified'")
    
    # 如果提供了override_eta_factor，则使用它
    if override_eta_factor is not None:
        eta_factor = override_eta_factor
    
    total_sigma_x_profile = np.zeros_like(y_coords_tunnel, dtype=float)

    for point in grouting_points_profile:
        grout_cylinder_height = point.z_bottom - point.z_top
        if grout_cylinder_height <= 0 or point.v_inj_point == 0:
            continue # Skip this point if height is invalid or no injection volume
        
        current_point_sigma_x_profile = np.zeros_like(y_coords_tunnel, dtype=float)
        x_src_grout = point.x_coord 
        y_src_grout = point.y_coord

        if point.non_uniform_expansion:
            # 使用 const.NON_UNIFORM_EXPANSION_FACTORS 和 const.NON_UNIFORM_SEGMENT_HEIGHT
            # 确保总高度与 grout_z_bottom - grout_z_top 匹配
            num_non_uniform_segments = len(const.NON_UNIFORM_EXPANSION_FACTORS)
            calculated_total_height = num_non_uniform_segments * const.NON_UNIFORM_SEGMENT_HEIGHT
            if abs(calculated_total_height - grout_cylinder_height) > 1e-3:
                print(f"Warning for grout point at x={point.x_coord},y={point.y_coord}: Non-uniform model height ({calculated_total_height}m) \
                      does not match grout interval height ({grout_cylinder_height}m). Using model height.")
            
            total_factor_sum = np.sum(const.NON_UNIFORM_EXPANSION_FACTORS)
            alpha_strength_total_for_V_inj = const.XI_INJ * point.v_inj_point # Use v_inj_point

            for i_segment in range(num_non_uniform_segments):
                # z_src_segment is the center of THIS non-uniform segment
                # Segments are from point.z_top downwards for this grout point
                z_src_segment_center = point.z_top + (i_segment + 0.5) * const.NON_UNIFORM_SEGMENT_HEIGHT
                segment_factor = const.NON_UNIFORM_EXPANSION_FACTORS[i_segment]
                alpha_strength_this_segment = (alpha_strength_total_for_V_inj * segment_factor) / total_factor_sum

                for i_obs in range(len(y_coords_tunnel)):
                    y_obs_tunnel = y_coords_tunnel[i_obs]
                    sigma_x_contrib = _calculate_sigma_x_from_point_source(
                        x_obs=tunnel_x_center,
                        y_obs=y_obs_tunnel,
                        z_obs=tunnel_z_center,
                        x_src=x_src_grout,
                        y_src=y_src_grout,
                        z_src=z_src_segment_center,
                        E=const.E_SOIL,
                        nu=const.NU_SOIL,
                        alpha_strength=alpha_strength_this_segment,
                        omit_complex_surface_terms=omit_complex_surface_terms
                    )
                    current_point_sigma_x_profile[i_obs] += sigma_x_contrib
        else: # Uniform expansion for this point
            num_segments_for_point = point.num_segments
            dz_segment = grout_cylinder_height / num_segments_for_point
            alpha_strength_per_segment = (const.XI_INJ * point.v_inj_point / grout_cylinder_height) * dz_segment
            
            for i_obs in range(len(y_coords_tunnel)):
                y_obs_tunnel = y_coords_tunnel[i_obs]
                current_sigma_at_y_for_point_segment = 0.0
                for k_segment in range(num_segments_for_point):
                    z_src_segment_center = point.z_top + (k_segment + 0.5) * dz_segment
                    sigma_x_contrib = _calculate_sigma_x_from_point_source(
                        x_obs=tunnel_x_center,
                        y_obs=y_obs_tunnel,
                        z_obs=tunnel_z_center,
                        x_src=x_src_grout,
                        y_src=y_src_grout,
                        z_src=z_src_segment_center,
                        E=const.E_SOIL,
                        nu=const.NU_SOIL,
                        alpha_strength=alpha_strength_per_segment,
                        omit_complex_surface_terms=omit_complex_surface_terms
                    )
                    current_sigma_at_y_for_point_segment += sigma_x_contrib
                current_point_sigma_x_profile[i_obs] = current_sigma_at_y_for_point_segment
        
        total_sigma_x_profile += current_point_sigma_x_profile

    return total_sigma_x_profile * eta_factor

if __name__ == '__main__':
    # --- 测试参数 (参照论文3.1节 工程案例) ---
    # 注浆参数 (使用新的 GroutPoint 结构)
    # Example 1: Single grout point, similar to previous default
    single_grout_point_profile = [const.DEFAULT_GROUT_POINT_FIG5B]
    # Update V_inj_example if needed, or use it from the point definition
    V_inj_example_total = sum(p.v_inj_point for p in single_grout_point_profile)

    # Example 2: Multiple grout points (Hypothetical)
    multiple_grout_points_profile = [
        GroutPoint(x_coord=const.R_TUNNEL + 10.3, y_coord=-2.0, z_top=15.0, z_bottom=20.0, v_inj_point=1.5, non_uniform_expansion=False, num_segments=30),
        GroutPoint(x_coord=const.R_TUNNEL + 10.3, y_coord=0.0,  z_top=15.0, z_bottom=20.0, v_inj_point=1.0, non_uniform_expansion=True), # num_segments ignored if non-uniform
        GroutPoint(x_coord=const.R_TUNNEL + 10.3, y_coord=2.0,  z_top=16.0, z_bottom=21.0, v_inj_point=1.5, non_uniform_expansion=False, num_segments=30),
    ]
    V_inj_example_total_multi = sum(p.v_inj_point for p in multiple_grout_points_profile)

    # 选择一个测试配置
    # current_grout_profile_for_test = single_grout_point_profile
    # V_inj_total_for_test = V_inj_example_total
    current_grout_profile_for_test = multiple_grout_points_profile
    V_inj_total_for_test = V_inj_example_total_multi

    # 隧道参数
    y_coords = const.Y_COORDS 
    z_tunnel = const.Z_TUNNEL_CENTER # 16.1 m
    x_tunnel = const.X_TUNNEL_CENTER # 0.0 m

    print(f"Calculating stress profile for V_inj = {V_inj_total_for_test} m^3")
    print(f"Tunnel X: {x_tunnel} m, Z: {z_tunnel} m")
    print(f"Y coordinates from {y_coords[0]} to {y_coords[-1]} m, {len(y_coords)} points")
    print(f"Soil E: {const.E_SOIL/1e6} MPa, nu: {const.NU_SOIL}")
    print(f"Eta correction (default): {const.ETA_CORRECTION}")
    print(f"Xi injection: {const.XI_INJ}")
    print(f"Testing with {len(current_grout_profile_for_test)} grout point(s), total V_inj = {V_inj_total_for_test} m^3")

    sigma_x_profile_example = calculate_total_sigma_x_profile(
        grouting_points_profile=current_grout_profile_for_test,
        y_coords_tunnel=y_coords,
        tunnel_x_center=const.X_TUNNEL_CENTER,
        tunnel_z_center=const.Z_TUNNEL_CENTER
        # omit_complex_surface_terms defaults to False
        # override_eta_factor defaults to None (uses const.ETA_CORRECTION)
    )

    print("\nCalculated sigma_x profile (Default Eta, Full Formula, Pa):")
    for i in range(len(y_coords)):
        if i % (len(y_coords)//10) == 0 or i == len(y_coords)-1: # Print a few values
             print(f"y = {y_coords[i]:.1f} m, sigma_x = {sigma_x_profile_example[i]:.2e} Pa")

    # Test non-uniform expansion - this is now controlled per GroutPoint
    # The main test above will use non-uniform if any point is set to it.
    # We can create a specific profile for testing non-uniform if needed.
    non_uniform_test_profile = [
        GroutPoint(const.R_TUNNEL + 10.3, 0.0, 15.0, 20.0, V_inj_total_for_test, non_uniform_expansion=True)
    ]
    if not any(p.non_uniform_expansion for p in current_grout_profile_for_test):
        print("\nCalculating stress profile (Explicit Non-Uniform Test, Default Eta, Full Formula)...")
        sigma_x_profile_non_uniform = calculate_total_sigma_x_profile(
            grouting_points_profile=non_uniform_test_profile,
            y_coords_tunnel=y_coords,
            tunnel_x_center=const.X_TUNNEL_CENTER,
            tunnel_z_center=const.Z_TUNNEL_CENTER
        )
        print("\nCalculated sigma_x profile (Non-Uniform Specific Test, Pa):")
        for i in range(len(y_coords)):
            if i % (len(y_coords)//10) == 0 or i == len(y_coords)-1:
                 print(f"y = {y_coords[i]:.1f} m, sigma_x = {sigma_x_profile_non_uniform[i]:.2e} Pa")

    # Test omitting terms with eta = 1.0
    print("\nCalculating stress profile (OMITTED TERMS, Eta=1.0)...")
    sigma_x_profile_omitted_eta1 = calculate_total_sigma_x_profile(
        grouting_points_profile=current_grout_profile_for_test, # Use the selected profile
        y_coords_tunnel=y_coords,
        tunnel_x_center=const.X_TUNNEL_CENTER,
        tunnel_z_center=const.Z_TUNNEL_CENTER,
        omit_complex_surface_terms=True,
        override_eta_factor=1.0
    )
    print("\nCalculated sigma_x profile (Omitted Terms, Eta=1.0, Pa):")
    for i in range(len(y_coords)):
        if i % (len(y_coords)//10) == 0 or i == len(y_coords)-1:
             print(f"y = {y_coords[i]:.1f} m, sigma_x = {sigma_x_profile_omitted_eta1[i]:.2e} Pa")

    # 尝试与论文中的图6或类似图的应力数量级进行粗略比较 (定性)
    # 论文中没有直接给出附加应力 σ'x 的分布图，只有位移图。
    # 不过，附加应力是中间步骤，其形态应与位移形态相关，中间大两边小。

    # 绘制示例图
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        plt.plot(y_coords, sigma_x_profile_example / 1000, label=f'Uniform V_inj = {V_inj_total_for_test} m^3')
        plt.plot(y_coords, sigma_x_profile_non_uniform / 1000, label=f'Non-Uniform V_inj = {V_inj_total_for_test} m^3', linestyle='--')
        plt.xlabel("y' (m)")
        plt.ylabel("附加水平应力 σ'x (kPa)")
        plt.title("计算得到的附加水平应力剖面 (修正简化法)")
        plt.grid(True)
        plt.legend()
        # plt.savefig("temp_stress_profile.png")
        plt.show()
    except ImportError:
        print("Matplotlib not installed. Cannot plot example stress profile.") 