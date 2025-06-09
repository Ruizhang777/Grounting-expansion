'''
plotting.py

包含生成和保存复现论文中各个图表的函数。
'''
import numpy as np
import matplotlib.pyplot as plt
import os
import constants as const # type: ignore

# Configure Matplotlib for Chinese characters
plt.rcParams['font.sans-serif'] = ['SimHei']  # Default to SimHei
plt.rcParams['axes.unicode_minus'] = False  # Ensure minus sign is displayed correctly

# Fallback fonts if SimHei is not available (Windows specific)
try:
    plt.figure() # Trigger font loading
    plt.close()
except Exception:
    try:
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
        plt.figure()
        plt.close()
        print("SimHei not found, using Microsoft YaHei.")
    except Exception:
        print("Neither SimHei nor Microsoft YaHei found. Chinese characters might not display correctly.")

def ensure_output_dir():
    if not os.path.exists(const.OUTPUT_DIR):
        os.makedirs(const.OUTPUT_DIR)
    print(f"Figures will be saved to: {os.path.abspath(const.OUTPUT_DIR)}")

def plot_fig6_and_similar(
    y_coords,
    calculated_displacement_mm, 
    main_calc_label="calculated distribution", # 改为英文标签
    main_calc_style='rs-', # 红色方形连线
    measured_data_mm=None, # Optional: for Fig 6, this is {y, displacement}
    finite_element_data_mm=None, # Optional: for Fig 6
    title="",  # 移除默认标题
    filename="fig_comparison.png",
    extra_lines=None # list of dicts: {'y': [], 'disp': [], 'label': str, 'style': str}
):
    ensure_output_dir()
    plt.figure(figsize=(10, 6.5))  # 恢复原来的图大小
    
    # 绘制计算分布 - 使用传入的main_calc_label和main_calc_style
    plt.plot(y_coords, calculated_displacement_mm, main_calc_style, label=main_calc_label, 
             linewidth=1.5, markersize=4)
    
    if measured_data_mm is not None and measured_data_mm.get('y_prime') is not None:
        # 实测数据使用蓝色菱形点匹配论文原图
        plt.plot(measured_data_mm['y_prime'], measured_data_mm['displacement'], 'D', 
                color='blue', label='measured data', markersize=6)
        
    if finite_element_data_mm is not None and finite_element_data_mm.get('y_prime') is not None:
        # 有限元数据使用绿色三角形标记连线匹配论文原图
        plt.plot(finite_element_data_mm['y_prime'], finite_element_data_mm['displacement'], '^-', 
                color='green', label='finite element data distribution', linewidth=1.5, markersize=5)

    if extra_lines:
        for line_data in extra_lines:
            style = line_data.get('style', '--')
            plt.plot(line_data['y'], line_data['disp'], style, label=line_data['label'], 
                    linewidth=1.5, markersize=4)
            
    plt.xlabel("y (m)", fontsize=12)  # 英文标签匹配论文原图
    plt.ylabel("Horizontal displacement (mm)", fontsize=12)  # 英文标签匹配论文原图
    # 移除标题以完全匹配论文原图
    # plt.title(title, fontsize=14)
    plt.xlim([const.Y_MIN, const.Y_MAX])
    plt.ylim([0, 4])  # 设置y轴范围为0-4mm匹配论文原图
    # 移除网格线以匹配论文原图
    # plt.grid(True, linestyle='--', alpha=0.7)
    
    # 调整图例位置和样式，匹配原论文图7
    plt.legend(fontsize=10, loc='upper right', frameon=False)
    
    plt.tight_layout() # Adjust layout to prevent clipping
    save_path = os.path.join(const.OUTPUT_DIR, filename)
    plt.savefig(save_path, dpi=const.FIGURE_DPI)
    print(f"Saved: {save_path}")
    plt.close()

def plot_fig9(
    y_coords,
    displacement_results, # dict: { "2 m^3": [], "4 m^3": [], ... }
    title="图9: 不同注浆体积对隧道水平位移补偿效果的影响",
    filename="fig9_influence_of_grouting_volume.png",
    legend_title="注浆体积"
):
    ensure_output_dir()
    plt.figure(figsize=(10, 6.5))  # 恢复原来的图大小
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    linestyles = ['-', '--', '-.', ':']
    
    # Sort keys by the numeric part for consistent plotting order
    try:
        # Attempt to sort numerically if keys follow "Value Unit" format
        sorted_keys = sorted(displacement_results.keys(), key=lambda x: float(str(x).split()[0]))
    except ValueError:
        # If conversion fails, use keys as they are (e.g., for Fig 14 style)
        sorted_keys = list(displacement_results.keys())

    for i, key in enumerate(sorted_keys):
        displacements = displacement_results[key]
        plt.plot(y_coords, displacements, 
                 color=colors[i % len(colors)], 
                 linestyle=linestyles[(i // len(colors)) % len(linestyles)], 
                 label=f"{key}", linewidth=1.5) # Changed label slightly for more general use
            
    plt.xlabel("y\' (m)", fontsize=12)
    plt.ylabel("Horizontal displacement (mm)", fontsize=12)
    plt.title(title, fontsize=14)
    plt.xlim([const.Y_MIN, const.Y_MAX])
    plt.ylim(bottom=0) # Based on Figure 9, displacement is positive and starts from 0
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10, title=legend_title)
    plt.tight_layout()
    save_path = os.path.join(const.OUTPUT_DIR, filename)
    plt.savefig(save_path, dpi=const.FIGURE_DPI)
    print(f"Saved: {save_path}")
    plt.close()

def plot_fig10(
    grouting_amounts, # array of V_inj values
    max_displacements, # array of corresponding max horizontal displacements
    theoretical_data=True, # True if it\'s the "Theoretical" line from Fig 10
    title="图10: 注浆体积与最大水平位移关系",
    filename="fig10_grouting_amount_vs_max_disp.png",
    xlabel_text="Grouting amount (m³)"
):
    ensure_output_dir()
    plt.figure(figsize=(8, 6))
    
    label = "Theoretical (修正简化法)" if theoretical_data else "Calculated Data"
    plt.plot(grouting_amounts, max_displacements, 'b-o', label=label, linewidth=1.5, markersize=5)
            
    plt.xlabel(xlabel_text, fontsize=12)
    plt.ylabel("Maximum horizontal displacement (mm)", fontsize=12)
    plt.title(title, fontsize=14)
    # plt.xlim(left=0) # Adjust as needed
    plt.ylim(bottom=0)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    save_path = os.path.join(const.OUTPUT_DIR, filename)
    plt.savefig(save_path, dpi=const.FIGURE_DPI)
    print(f"Saved: {save_path}")
    plt.close()

# Similar functions for Fig11, Fig12, Fig14 will be needed.
# plot_fig11: Influence of grouting distance (similar to plot_fig9)
# plot_fig12: Relationship between grouting distance and max displacement (similar to plot_fig10)
# plot_fig14: Influence of grouting depth (similar to plot_fig9, but x-axis is y', legend for exp groups)

if __name__ == '__main__':
    print("Testing plotting.py")
    ensure_output_dir()

    # Test plot_fig6_and_similar (like Fig 7)
    y_coords_test = const.Y_COORDS
    calc_disp_test = 3.5 * np.exp(-(y_coords_test / 10)**2) + 0.1 * np.sin(y_coords_test/5)
    simplified_disp_test = 3.0 * np.exp(-(y_coords_test / 10.5)**2)
    
    extra_lines_fig7 = [
        {'y': y_coords_test, 'disp': simplified_disp_test, 'label': 'Simplified method (example)', 'style': 'g--'}
    ]
    plot_fig6_and_similar(y_coords_test, calc_disp_test, 
                            main_calc_label="主计算线 (示例)",
                            main_calc_style='r--', # Test the new style parameter
                            measured_data_mm=const.MEASURED_DATA_FIG7,
                            title="测试图7: 水平位移比较 (示例数据)",
                            filename="test_fig7_example.png",
                            extra_lines=extra_lines_fig7)

    # Test plot_fig9
    disp_results_test_fig9 = {
        "2 m$^3$": 2.0 * np.exp(-(y_coords_test / 12)**2),
        "4 m$^3$": 4.0 * np.exp(-(y_coords_test / 10)**2),
        "6 m$^3$": 6.0 * np.exp(-(y_coords_test / 9)**2),
        "8 m$^3$": 7.5 * np.exp(-(y_coords_test / 8.5)**2),
        "10 m$^3$": 8.0 * np.exp(-(y_coords_test / 8)**2),
    }
    plot_fig9(y_coords_test, disp_results_test_fig9, filename="test_fig9_example.png")

    # Test plot_fig10
    grouting_amounts_test = np.array([2, 4, 6, 8, 10])
    max_disps_test = np.array([v*0.8 for v in grouting_amounts_test]) # Example data
    plot_fig10(grouting_amounts_test, max_disps_test, filename="test_fig10_example.png")

    print("Plotting tests complete. Check the 'figures' directory.") 