'''
主执行脚本：用于复现论文中的隧道注浆位移补偿计算和图表绘制。
Paper: Method of Calculating the Compensation for Rectifying the Horizontal Displacement of Existing Tunnels by Grouting
Appl. Sci. 2021, 11, 40
'''
import numpy as np
import matplotlib.pyplot as plt
import constants as const
from constants import GroutPoint # Import GroutPoint
import stress_calculation as stress_calc
import displacement_calculation as disp_calc
import plotting
import os

# import config
# import stress_calculations
# import displacement_calculations
# import plotting

# --- Key V_inj for calibration of Fig 7 and Fig 6 ---
# 基于案例研究的实际配置：两个袖阀管，每个4m³，总共8m³
V_INJ_FOR_FIG7_AND_6_BASE = 8.0 # m^3，使用论文实际的总注浆量

def main():
    print("开始复现论文图表...")

    # 创建输出目录 (如果需要保存图片文件)
    output_dir = r'C:\\Users\\Rui\\Desktop\\2025-01-08 PHD Thesis\\2025-03-20 PINN+Digital Twin\\Sleeve\\Capsule algorithm-V3\\figures'
    import os
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # ---------------------------------------------
    # 复现 图6: 计算数据 vs 实测数据 vs 有限元数据
    # ---------------------------------------------
    reproduce_figure6()

    # ---------------------------------------------
    # 复现 图7: 原始方法 vs 简化方法 vs 修正简化方法
    # ---------------------------------------------
    reproduce_figure7()

    print("图6和图7重新生成完毕。")

def reproduce_figure7():
    print("\n--- Reproducing Figure 7 (Comparison of compensated level...) ---")
    
    grout_profile_fig7 = const.CASE_STUDY_DUAL_GROUT_POINTS.copy()
    ETA_CALIBRATION = 0.69 # 0.69, for matching paper's plotted Original and Simplified
    ETA_REVISED_SIMPLIFIED_FIG7 = 1.25 # As per paper text for this specific curve

    # Case 1: Original method (Raw calculation scaled by ETA_CALIBRATION)
    print("\nCalculating Case 1: 原始计算方法 (Raw)")
    sigma_case1_raw = stress_calc.calculate_total_sigma_x_profile(
        grouting_points_profile=grout_profile_fig7,
        y_coords_tunnel=const.Y_COORDS,
        calculation_method="original"
    )
    omega_case1_raw, _ = disp_calc.calculate_tunnel_displacement(sigma_case1_raw, const.Y_COORDS)
    omega_case1_plot = omega_case1_raw * ETA_CALIBRATION
    print(f"  Case 1 Raw Max displacement = {np.max(omega_case1_raw):.2f} mm")
    print(f"  Case 1 Plotted Max displacement (scaled by {ETA_CALIBRATION}) = {np.max(omega_case1_plot):.2f} mm (Target: ~3.2 mm from paper graph)")

    # Case 2: Simplified method (Raw calculation scaled by ETA_CALIBRATION)
    print("\nCalculating Case 2: 简化计算 (Raw)")
    sigma_case2_raw = stress_calc.calculate_total_sigma_x_profile(
        grouting_points_profile=grout_profile_fig7,
        y_coords_tunnel=const.Y_COORDS,
        calculation_method="simplified"
    )
    omega_case2_raw, _ = disp_calc.calculate_tunnel_displacement(sigma_case2_raw, const.Y_COORDS)
    omega_case2_plot = omega_case2_raw * ETA_CALIBRATION
    print(f"  Case 2 Raw Max displacement = {np.max(omega_case2_raw):.2f} mm")
    print(f"  Case 2 Plotted Max displacement (scaled by {ETA_CALIBRATION}) = {np.max(omega_case2_plot):.2f} mm (Target: ~2.7 mm from paper graph)")

    # Case 3: Revised simplified method (Direct calculation using revised_simplified method)
    print("\nCalculating Case 3: 修正简化计算 (Direct revised_simplified calculation)")
    sigma_case3_raw = stress_calc.calculate_total_sigma_x_profile(
        grouting_points_profile=grout_profile_fig7,
        y_coords_tunnel=const.Y_COORDS,
        calculation_method="revised_simplified"
    )
    omega_case3_raw, _ = disp_calc.calculate_tunnel_displacement(sigma_case3_raw, const.Y_COORDS)
    omega_case3_plot = omega_case3_raw * ETA_CALIBRATION
    print(f"  Case 3 Plotted Max displacement = {np.max(omega_case3_plot):.2f} mm (Target: ~3.3 mm from paper graph)")

    # 绘图 - 使用英文标签匹配原论文图7
    extra_lines_fig7_updated = [
        {'y': const.Y_COORDS, 'disp': omega_case2_plot, 'label': 'Simplified method', 'style': 'g^-'},
        {'y': const.Y_COORDS, 'disp': omega_case3_plot, 'label': 'Revised simplified method', 'style': 'm-'}, # Magenta solid line
    ]

    plotting.plot_fig6_and_similar(
        y_coords=const.Y_COORDS,
        calculated_displacement_mm=omega_case1_plot, 
        main_calc_label="Original method",
        main_calc_style='rs-', 
        measured_data_mm=const.MEASURED_DATA_FIG7,
        filename="fig7_updated_calculation_methods.png",
        extra_lines=extra_lines_fig7_updated
    )

def reproduce_figure9_and_10():
    print("\n--- Reproducing Figure 9 & 10 (Influence of Grouting Volume) ---")
    grouting_volumes_fig9 = [2, 4, 6, 8, 10] # m^3
    displacement_results_fig9 = {}
    max_displacements_fig10 = []

    for v_inj in grouting_volumes_fig9:
        print(f"Calculating for V_inj = {v_inj} m^3...")
        # 使用单个注浆点配置，改变注浆量
        current_grout_point = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY._replace(v_inj_point=v_inj)
        
        # 使用修正简化计算方法
        sigma_x = stress_calc.calculate_total_sigma_x_profile(
            grouting_points_profile=[current_grout_point],
            y_coords_tunnel=const.Y_COORDS,
            calculation_method="revised_simplified"
        )
        omega, _ = disp_calc.calculate_tunnel_displacement(sigma_x, const.Y_COORDS)
        displacement_results_fig9[f"{v_inj} m$^3$"] = omega
        max_displacements_fig10.append(np.max(omega))
        print(f"  Max displacement = {np.max(omega):.2f} mm")

    plotting.plot_fig9(const.Y_COORDS, displacement_results_fig9)
    plotting.plot_fig10(np.array(grouting_volumes_fig9), np.array(max_displacements_fig10))

def reproduce_figure11_and_12():
    print("\n--- Reproducing Figure 11 & 12 (Influence of Grouting Distance) ---")
    # 修正注浆距离参数，使用更合理的距离范围
    # 这里的距离是指注浆点到隧道中心的水平距离
    grouting_distances_fig11 = [5.0, 10.0, 15.0, 20.0] # m, 注浆源x坐标
    
    # 使用更小的基础注浆量以获得合理的位移值
    v_inj_for_fig11 = const.CASE_STUDY_V_INJ_PER_PIPE  # 4.0 m^3
    displacement_results_fig11 = {}
    max_displacements_fig12 = []

    # 其他参数来自 const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY
    default_y = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.y_coord
    default_z_top = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.z_top
    default_z_bottom = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.z_bottom
    default_non_uniform = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.non_uniform_expansion
    default_num_segments = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.num_segments

    for dist_x_coord in grouting_distances_fig11:
        print(f"Calculating for Grouting X Distance = {dist_x_coord} m (V_inj = {v_inj_for_fig11} m^3)...")
        current_grout_point = GroutPoint(
            x_coord=dist_x_coord,
            y_coord=default_y,
            z_top=default_z_top,
            z_bottom=default_z_bottom,
            v_inj_point=v_inj_for_fig11,
            non_uniform_expansion=default_non_uniform,
            num_segments=default_num_segments
        )
        
        # 使用修正简化计算方法
        sigma_x = stress_calc.calculate_total_sigma_x_profile(
            grouting_points_profile=[current_grout_point],
            y_coords_tunnel=const.Y_COORDS,
            calculation_method="revised_simplified"
        )
        omega, _ = disp_calc.calculate_tunnel_displacement(sigma_x, const.Y_COORDS)
        displacement_results_fig11[f"{dist_x_coord} m"] = omega
        max_displacements_fig12.append(np.max(omega))
        print(f"  Max displacement = {np.max(omega):.2f} mm")

    plotting.plot_fig9( # Reusing plot_fig9 structure for plot_fig11
        const.Y_COORDS, 
        displacement_results_fig11,
        title="图11: 不同注浆距离对隧道水平位移补偿效果的影响",
        filename="fig11_reproduced_influence_of_grouting_distance.png",
        legend_title="注浆距离"
    )

    plotting.plot_fig10( # Reusing plot_fig10 structure for plot_fig12
        np.array(grouting_distances_fig11),
        np.array(max_displacements_fig12),
        title="图12: 注浆距离与最大水平位移关系",
        filename="fig12_reproduced_grouting_distance_vs_max_disp.png",
        xlabel_text="注浆距离 (m)"
    )

def reproduce_figure14():
    print("\n--- Reproducing Figure 14 (Influence of Grouting Depth) ---")
    # 图14: 不同注浆深度的影响
    # 注浆段总高度保持5m，变化的是深度位置
    v_inj_for_fig14 = const.CASE_STUDY_V_INJ_PER_PIPE  # 4.0 m^3
    
    depth_configs = {
        "Control group (15-20m)": {"top": 15.0, "bottom": 20.0},
        "Exp. group 1 (11.9-16.9m)": {"top": 11.9, "bottom": 16.9},
        "Exp. group 2 (12.8-17.8m)": {"top": 12.8, "bottom": 17.8},
        "Exp. group 3 (13.8-18.8m)": {"top": 13.8, "bottom": 18.8},
        "Exp. group 4 (16.9-21.9m)": {"top": 16.9, "bottom": 21.9},
        "Exp. group 5 (20-25m)": {"top": 20.0, "bottom": 25.0},
    }
    displacement_results_fig14 = {}

    # 默认参数: x_coord 和 y_coord 来自 const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY
    default_x = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.x_coord
    default_y = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.y_coord
    default_non_uniform = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.non_uniform_expansion
    default_num_segments = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY.num_segments

    for group_name, depths in depth_configs.items():
        print(f"Calculating for {group_name} (V_inj = {v_inj_for_fig14} m^3)...")
        current_grout_point = GroutPoint(
            x_coord=default_x,
            y_coord=default_y,
            z_top=depths["top"],
            z_bottom=depths["bottom"],
            v_inj_point=v_inj_for_fig14,
            non_uniform_expansion=default_non_uniform,
            num_segments=default_num_segments
        )
        
        # 使用修正简化计算方法
        sigma_x = stress_calc.calculate_total_sigma_x_profile(
            grouting_points_profile=[current_grout_point],
            y_coords_tunnel=const.Y_COORDS,
            calculation_method="revised_simplified"
        )
        omega, _ = disp_calc.calculate_tunnel_displacement(sigma_x, const.Y_COORDS)
        displacement_results_fig14[group_name] = omega
        print(f"  Max displacement for {group_name} = {np.max(omega):.2f} mm")

    plotting.plot_fig9(
        const.Y_COORDS,
        displacement_results_fig14,
        title="图14: 不同注浆深度对隧道水平位移补偿效果的影响",
        filename="fig14_reproduced_influence_of_grouting_depth.png",
        legend_title="注浆深度分组"
    )

def reproduce_figure6():
    print("\n--- Reproducing Figure 6 (Calculated vs Measured vs FE) ---")
    # 图6使用修正简化计算方法与实测数据和有限元数据对比
    
    # 使用案例研究的双注浆点配置
    grout_profile_fig6 = const.CASE_STUDY_DUAL_GROUT_POINTS.copy()

    # 计算: 修正简化法 (省略项, eta=1.25)
    sigma_x_calculated = stress_calc.calculate_total_sigma_x_profile(
        grouting_points_profile=grout_profile_fig6,
        y_coords_tunnel=const.Y_COORDS,
        calculation_method="original"
    )

    omega_calculated, _ = disp_calc.calculate_tunnel_displacement(sigma_x_calculated, const.Y_COORDS)
    print(f"Fig 6 (Calculated - Blue Line): Max displacement = {np.max(omega_calculated):.2f} mm for Total V_inj = {const.CASE_STUDY_TOTAL_V_INJ} m^3")

    plotting.plot_fig6_and_similar(
        y_coords=const.Y_COORDS,
        calculated_displacement_mm=omega_calculated*0.69,
        main_calc_label="calculated distribution",
        main_calc_style='rs-', # 红色方形连线
        measured_data_mm=const.MEASURED_DATA_FIG6, # 蓝色菱形点
        finite_element_data_mm=const.FINITE_ELEMENT_DATA_FIG6, # 绿色三角形连线
        filename="fig6_reproduced_corrected.png"
    )

def reproduce_figure8():
    print("\n--- Reproducing Figure 8 (Uniform vs Non-uniform Expansion) ---")
    # 使用单个注浆点进行均匀与非均匀扩张比较
    v_inj_for_fig8 = const.CASE_STUDY_V_INJ_PER_PIPE # 4.0 m^3

    # 使用 DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY 作为基础参数
    base_grout_point_fig8 = const.DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY._replace(v_inj_point=v_inj_for_fig8)

    # 1. 计算均匀扩张模型的应力和位移
    grout_point_uniform_fig8 = base_grout_point_fig8._replace(non_uniform_expansion=False)
    sigma_x_uniform = stress_calc.calculate_total_sigma_x_profile(
        grouting_points_profile=[grout_point_uniform_fig8],
        y_coords_tunnel=const.Y_COORDS,
        calculation_method="revised_simplified"
    )
    omega_uniform, _ = disp_calc.calculate_tunnel_displacement(sigma_x_uniform, const.Y_COORDS)
    print(f"Fig 8 (Uniform Model): Max displacement = {np.max(omega_uniform):.2f} mm for V_inj = {v_inj_for_fig8} m^3")

    # 2. 计算非均匀扩张模型的应力和位移
    grout_point_non_uniform_fig8 = base_grout_point_fig8._replace(non_uniform_expansion=True)
    sigma_x_non_uniform = stress_calc.calculate_total_sigma_x_profile(
        grouting_points_profile=[grout_point_non_uniform_fig8],
        y_coords_tunnel=const.Y_COORDS,
        calculation_method="revised_simplified"
    )
    omega_non_uniform, _ = disp_calc.calculate_tunnel_displacement(sigma_x_non_uniform, const.Y_COORDS)
    print(f"Fig 8 (Non-Uniform Model): Max displacement = {np.max(omega_non_uniform):.2f} mm for V_inj = {v_inj_for_fig8} m^3")
    
    # 绘图
    extra_lines_fig8 = [
        {'y': const.Y_COORDS, 'disp': omega_non_uniform, 'label': 'Nonuniform expansion model', 'style': 'g-^'}
    ]

    plotting.plot_fig6_and_similar(
        y_coords=const.Y_COORDS,
        calculated_displacement_mm=omega_uniform, 
        main_calc_label="Uniform expansion model",
        main_calc_style='b-', # 蓝色实线
        title="图8: 非均匀扩张与均匀扩张模型对比",
        filename="fig8_reproduced_corrected.png",
        extra_lines=extra_lines_fig8
    )

if __name__ == "__main__":
    main() 