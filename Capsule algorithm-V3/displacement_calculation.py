'''
displacement_calculation.py

包含基于能量法和傅里叶级数计算隧道水平位移的函数 (论文式21-26)。
'''
import numpy as np
import constants as const # type: ignore
from scipy.integrate import simps # Simpson's rule for numerical integration

def calculate_fourier_coefficients_of_stress(
    sigma_x_profile, 
    y_coords # Should be from Y_MIN to Y_MAX
):
    '''
    计算附加应力剖面 sigma_x_profile 的傅里叶余弦级数系数。
    σ_x(y) ≈ A0/2 + Σ Am * cos(m * π * (y - Y_MIN) / L) for m=1...N
    Or, σ_x(y) ≈ Σ sig_coeffs[m] * cos(m * π * y_shifted / L) for m=0...N
    where y_shifted = y - Y_MIN, L = ANALYSIS_LENGTH_L

    返回: {σx}T 向量 (论文中的 notation)
    '''
    L = const.ANALYSIS_LENGTH_L
    N_terms = const.N_FOURIER_TERMS
    y_shifted = y_coords - const.Y_MIN # y' from 0 to L

    sigma_coeffs = np.zeros(N_terms)

    # Zeroth coefficient (A0 or sigma_coeffs[0])
    # A0 = (1/L) * integral(sigma_x(y') dy') from 0 to L for basis {1}
    # If basis is {1, cos(πy'/L), ...}, then T0(y')=1. Integral T0^2 = L
    # sig_coeffs[0] = (1/L) * simps(sigma_x_profile, y_shifted)
    # The paper formulation AT = ([Kc] + [Ks])^-1 {σx}T uses {σx}T where Tn(y) = {1, cos(πy/L), ...}
    # Let's assume {σx}T_m = integral( σ_x(y) * Tm(y) dy ) / (D_tunnel) ? No, this is force resultant.
    # The paper's energy method: dW_ext / da_m = {σx}T_m
    # W_ext = integral (σx(y) * D_tunnel * w(y)) dy (assuming σx is stress, D_tunnel*w(y) is area deformation)
    # w(y) = sum a_m Tm(y). So dW_ext/da_m = integral( σx(y) * D_tunnel * Tm(y) ) dy
    # So, sigma_coeffs[m] = D_tunnel * integral(sigma_x_profile(y) * cos(m*π*y'/L) dy ) from y'=0 to L
    # (where m=0 term is just integral(sigma_x_profile(y) * D_tunnel dy))
    
    # For T0(y') = 1
    integrand0 = sigma_x_profile * const.D_TUNNEL
    sigma_coeffs[0] = simps(integrand0, y_shifted)
    
    # For Tm(y') = cos(m*π*y'/L), m > 0
    for m in range(1, N_terms):
        basis_func_m = np.cos(m * np.pi * y_shifted / L)
        integrand_m = sigma_x_profile * const.D_TUNNEL * basis_func_m
        sigma_coeffs[m] = simps(integrand_m, y_shifted)
        
    return sigma_coeffs

def calculate_tunnel_displacement(
    sigma_x_profile, 
    y_coords, # Should be from Y_MIN to Y_MAX
    include_shear_rigidity: bool = True  # 新增：是否包含环间剪切刚度Kr项
):
    '''
    根据给定的附加应力剖面计算隧道水平位移。

    参数:
        sigma_x_profile: 附加水平应力剖面 (Pa) (array)
        y_coords: 对应的y坐标 (m) (array)
        include_shear_rigidity: 是否包含环间剪切刚度Kr项

    返回:
        omega_y: 隧道水平位移 (mm) (array)
        y_coords: 输入的y坐标 (for plotting convenience)
    '''
    N_terms = const.N_FOURIER_TERMS
    L = const.ANALYSIS_LENGTH_L
    y_shifted = y_coords - const.Y_MIN # y' from 0 to L

    # 1. 计算应力傅里叶系数 {σx}T
    sigma_fourier_coeffs = calculate_fourier_coefficients_of_stress(sigma_x_profile, y_coords)

    # 2. 构建刚度矩阵 [Kc] 和 [Ks] 和 [Kr] (对角矩阵)
    Kc_diag = np.zeros(N_terms)  # 隧道抗弯刚度
    Ks_diag = np.zeros(N_terms)  # 地基反力刚度
    Kr_diag = np.zeros(N_terms)  # 环间剪切刚度 (新增)

    # Kc_m: 来自隧道抗弯刚度 EI_TUNNEL
    # d^2(cos(mπy'/L))/dy'^2 = -(mπ/L)^2 cos(mπy'/L)
    # Integral_0^L [ (mπ/L)^2 cos(mπy'/L) ]^2 dy' = (mπ/L)^4 * (L/2) for m > 0
    # Integral_0^L [ (0) ]^2 dy' = 0 for m = 0 (constant displacement term)
    # So Kc_0 = 0 (no resistance to uniform translation from bending stiffness)
    Kc_diag[0] = 1e-12 # Small epsilon to avoid singularity if other stiffness terms are also zero
    for m in range(1, N_terms):
        Kc_diag[m] = const.EI_TUNNEL * (m * np.pi / L)**4 * (L / 2.0)
    
    # Kr_m: 来自环间剪切刚度 KC_TUNNEL_SHEAR_RIGIDITY_GA (新增Kr项)
    # 基于梁的剪切变形理论，剪切刚度矩阵项为 GA * (m*π/L)^2
    # Integral_0^L [ (mπ/L) cos(mπy'/L) ]^2 dy' = (mπ/L)^2 * (L/2) for m > 0
    if include_shear_rigidity and const.KC_TUNNEL_SHEAR_RIGIDITY_GA > 1e-9:
        Kr_diag[0] = 1e-12  # 常位移项无剪切阻抗
        for m in range(1, N_terms):
            Kr_diag[m] = const.KC_TUNNEL_SHEAR_RIGIDITY_GA * (m * np.pi / L)**2 * (L / 2.0)
    
    # Ks_m: 来自地基反力 K_SUBGRADE_REACTION (Pa/m, force per area per displacement)
    # Effective spring constant per unit length of tunnel = K_SUBGRADE_REACTION * D_TUNNEL
    # Integral_0^L [ cos(mπy'/L) ]^2 dy' = L/2 for m > 0, and L for m = 0 (basis func = 1)
    ks_effective_per_length = const.K_SUBGRADE_REACTION * const.D_TUNNEL
    
    Ks_diag[0] = ks_effective_per_length * L # For T0=1, integral(T0^2) = L
    for m in range(1, N_terms):
        Ks_diag[m] = ks_effective_per_length * (L / 2.0)

    # 总刚度矩阵：[K_total] = [Kc] + [Ks] + [Kr]
    K_total_diag = Kc_diag + Ks_diag + Kr_diag
    
    # Handle potential division by zero if a K_total_diag term is zero, though unlikely with Ks.
    K_total_inv_diag = np.zeros_like(K_total_diag)
    for m in range(N_terms):
        if abs(K_total_diag[m]) < 1e-9:
            K_total_inv_diag[m] = 0 # Or a very large number if corresponding sigma_coeff is non-zero (problematic)
            if abs(sigma_fourier_coeffs[m]) > 1e-9:
                print(f"Warning: Near-zero stiffness K_total_diag[{m}] with non-zero force sigma_fourier_coeffs[{m}]")
        else:
            K_total_inv_diag[m] = 1.0 / K_total_diag[m]

    # 3. 求解位移的傅里叶系数 a_m (论文中的 AT)
    # AT = ([Kc] + [Ks] + [Kr])^-1 * {σx}T
    displacement_fourier_coeffs = K_total_inv_diag * sigma_fourier_coeffs

    # 4. 合成隧道位移 ω(y)
    # ω(y) = Σ a_m * Tm(y)
    omega_y = np.zeros_like(y_coords, dtype=float)
    
    # Add m=0 term (uniform displacement)
    omega_y += displacement_fourier_coeffs[0] * 1.0 # T0(y') = 1
    
    # Add m>0 terms
    for m in range(1, N_terms):
        basis_func_m = np.cos(m * np.pi * y_shifted / L)
        omega_y += displacement_fourier_coeffs[m] * basis_func_m

    return omega_y * 1000.0, y_coords # Convert to mm


if __name__ == '__main__':
    # --- 测试参数 ---
    print("Testing displacement_calculation.py")
    
    # 构造一个简单的对称应力剖面 (e.g., Gaussian-like for testing)
    y_test = const.Y_COORDS
    sigma_max_kPa = 50 # kPa
    sigma_test_profile = (sigma_max_kPa * 1000) * np.exp(- (y_test / (const.ANALYSIS_LENGTH_L/6))**2) 
    
    print(f"Test y_coords: {y_test.shape}, from {y_test[0]} to {y_test[-1]}")
    print(f"Test sigma_profile (max): {np.max(sigma_test_profile)/1000:.2f} kPa")
    print(f"EI_TUNNEL: {const.EI_TUNNEL:.2e} Nm^2")
    print(f"K_SUBGRADE_REACTION: {const.K_SUBGRADE_REACTION:.2e} Pa/m")
    print(f"KC_TUNNEL_SHEAR_RIGIDITY_GA: {const.KC_TUNNEL_SHEAR_RIGIDITY_GA:.2e} N")
    print(f"D_TUNNEL: {const.D_TUNNEL} m")
    print(f"N_FOURIER_TERMS: {const.N_FOURIER_TERMS}")
    print(f"ANALYSIS_LENGTH_L: {const.ANALYSIS_LENGTH_L} m")

    omega_profile, _ = calculate_tunnel_displacement(sigma_test_profile, y_test)

    print("\nCalculated omega_y profile (mm):")
    for i in range(len(y_test)):
        if i % (len(y_test)//10) == 0 or i == len(y_test)-1:
            print(f"y = {y_test[i]:.1f} m, omega = {omega_profile[i]:.4f} mm")

    # 绘图
    try:
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(10, 6))

        color = 'tab:red'
        ax1.set_xlabel("y' (m)")
        ax1.set_ylabel("附加水平应力 σ'x (kPa)", color=color)
        ax1.plot(y_test, sigma_test_profile / 1000, color=color, linestyle='--')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, linestyle=':')

        ax2 = ax1.twinx() # 共享x轴
        color = 'tab:blue'
        ax2.set_ylabel("隧道水平位移 ω(y) (mm)", color=color)
        ax2.plot(y_test, omega_profile, color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        plt.title("测试: 隧道水平位移计算")
        fig.tight_layout() # 调整布局以防止标签重叠
        # plt.savefig("temp_displacement_calc_test.png")
        plt.show()
    except ImportError:
        print("Matplotlib not installed. Cannot plot example displacement profile.") 