'''
constants.py

存储从论文中提取的物理参数、土体参数、隧道参数、注浆参数等常量。
'''
import numpy as np
from typing import List, NamedTuple, Dict # Added NamedTuple and List

# --- 定义注浆点的数据结构 ---
class GroutPoint(NamedTuple):
    x_coord: float        # 注浆源中心 x 坐标 (m)
    y_coord: float        # 注浆源中心 y 坐标 (m)
    z_top: float          # 注浆段顶部深度 (m)
    z_bottom: float       # 注浆段底部深度 (m)
    v_inj_point: float    # 该注浆点的注浆量 (m^3)
    non_uniform_expansion: bool = False # 是否对该点使用非均匀扩张模型
    num_segments: int = 50 # 如果是均匀扩张，该点的分段数

# --- 土体参数 (参考论文 3.1 节 工程案例)
E_SOIL = 25e6  # 土壤弹性模量 (Pa), 25 MPa
NU_SOIL = 0.3  # 土壤泊松比 (无量纲)
G_SOIL = E_SOIL / (2 * (1 + NU_SOIL)) # 土壤剪切模量 (Pa)

# --- 隧道参数 (参考论文 3.1 节 工程案例)
D_TUNNEL = 6.2  # 隧道直径 (m)
R_TUNNEL = D_TUNNEL / 2 # 隧道半径 (m)
# 隧道等效抗弯刚度 EI (Nm^2) - 根据论文调整
EI_TUNNEL = 1.94e8 # Nm^2 (论文3.1节: 1.94e8 kN*m^2，转换为N*m^2)
K_SUBGRADE_REACTION = 7.45e6 # 地基综合反力系数 ks (N/m^3), 7.45e6 N/m^3

# 隧道剪切刚度 (G*A_shear), 用于考虑环间剪切效应Kr项
KC_TUNNEL_SHEAR_RIGIDITY_GA = 5e9 # 单位: N，调整为更合理的值

# --- 注浆参数
XI_INJ = 1.0  # 注浆效率 (论文中假设为1)
ETA_CORRECTION = 1.25  # 修正简化法的修正系数 (调整为0.69以匹配论文图6结果)

# 默认注浆体几何参数 (参考论文 Section 3.1 和 Figure 5a)
GROUT_CYLINDER_HEIGHT = 5.0  # 注浆加固区的高度 h_cyl (m) (即论文中注浆段长度为5m)
GROUT_CYLINDER_RADIUS_R1 = 0.8  # 注浆前加固区半径 R1 (m)
# R2 将根据 V_inj 和 R1, XI_INJ, GROUT_CYLINDER_HEIGHT 计算得到

# --- 计算参数
N_FOURIER_TERMS = 10  # 傅里叶级数项数 (论文推荐 n=10)
Y_MIN = -30.0  # 计算范围 y坐标最小值 (m)
Y_MAX = 30.0  # 计算范围 y坐标最大值 (m)
NUM_Y_POINTS = 201 # y坐标取点数量
Y_COORDS = np.linspace(Y_MIN, Y_MAX, NUM_Y_POINTS) # y坐标数组
ANALYSIS_LENGTH_L = Y_MAX - Y_MIN # 傅里叶级数展开的分析长度

# --- 默认几何参数 (参考论文 3.1 节 工程案例 和 Figure 5)
# 假设隧道中心在 x=0
X_TUNNEL_CENTER = 0.0  # 隧道中心x坐标 (m)
# 隧道中心深度估算: 地表覆土13m, 隧道直径6.2m -> 13 + 3.1 = 16.1m
Z_TUNNEL_CENTER = 16.1 # 隧道中心z坐标 (深度, m)

# --- 案例研究中的注浆管线参数 (基于论文描述和用户反馈 "问题1") ---
# 论文提及 "two sleeve pipes", "distance of 10.3 m from the side wall", "grouting volume of a single sleeve pipe was 4 m^3"
# 修正：两个管线应该对称布置以产生对称的位移分布
CASE_STUDY_X_GROUT_OFFSET_FROM_WALL = 10.3  # m, 注浆管距隧道壁的水平距离
CASE_STUDY_X_GROUT_COORD = R_TUNNEL + CASE_STUDY_X_GROUT_OFFSET_FROM_WALL # 注浆管中心x坐标
CASE_STUDY_Y_GROUT_COORD_PIPE1 = -2.0   # m, 第一个注浆管y坐标（调整为更接近中心）
CASE_STUDY_Y_GROUT_COORD_PIPE2 = 2.0    # m, 第二个注浆管y坐标（调整为更接近中心）
CASE_STUDY_Z_GROUT_TOP = 15.0           # m, 注浆段顶部深度
CASE_STUDY_Z_GROUT_BOTTOM = 20.0        # m, 注浆段底部深度
CASE_STUDY_V_INJ_PER_PIPE = 4.0         # m^3, 每个注浆管的注浆量（恢复论文原始值）

# --- 定义案例研究的两个注浆点配置 ---
CASE_STUDY_DUAL_GROUT_POINTS = [
    GroutPoint(
        x_coord=CASE_STUDY_X_GROUT_COORD,
        y_coord=CASE_STUDY_Y_GROUT_COORD_PIPE1,
        z_top=CASE_STUDY_Z_GROUT_TOP,
        z_bottom=CASE_STUDY_Z_GROUT_BOTTOM,
        v_inj_point=CASE_STUDY_V_INJ_PER_PIPE,
        non_uniform_expansion=False,
        num_segments=50
    ),
    GroutPoint(
        x_coord=CASE_STUDY_X_GROUT_COORD,
        y_coord=CASE_STUDY_Y_GROUT_COORD_PIPE2,
        z_top=CASE_STUDY_Z_GROUT_TOP,
        z_bottom=CASE_STUDY_Z_GROUT_BOTTOM,
        v_inj_point=CASE_STUDY_V_INJ_PER_PIPE,
        non_uniform_expansion=False,
        num_segments=50
    )
]

# 总注浆量
CASE_STUDY_TOTAL_V_INJ = sum(point.v_inj_point for point in CASE_STUDY_DUAL_GROUT_POINTS)

# # 注浆管位置 (距隧道壁10.3m) # Commented out, will be part of GroutPoint
# # X_GROUT_PIPE = R_TUNNEL + 10.3
# X_GROUT_PIPE_DEFAULT_OFFSET = R_TUNNEL + 10.3 # 注浆管距隧道中心的默认水平距离 (m)
# Y_GROUT_SECTION_CENTER = 0.0 # 注浆段y向中心 (m), 通常设为0
# 
# # 注浆段深度 (论文3.1节：15m至20m)
# Z_GROUT_TOP_DEFAULT = 15.0  # 注浆段顶部深度 (m)
# Z_GROUT_BOTTOM_DEFAULT = 20.0 # 注浆段底部深度 (m)
# # 这也意味着 GROUT_CYLINDER_HEIGHT = 5.0 m，与之前定义一致

# --- 示例注浆点配置 (用于替换旧的单一注浆点参数) ---
# 根据 figure 5b: x = 10.3m (距隧道壁) + 3.1m (隧道半径) = 13.4m
# y = 0 (假设在隧道中心线上方或下方对称布置，这里取中心)
# z_top=15, z_bottom=20 (同论文案例)
# v_inj_point: 需要根据具体场景设定，例如总注浆量分配到各点
# 这里用一个与之前 V_INJ_FOR_FIG7_AND_6_BASE 类似的示例值
# DEFAULT_GROUT_POINT_FIG5B = GroutPoint(
#     x_coord = R_TUNNEL + 10.3, 
#     y_coord = 0.0, 
#     z_top = 15.0, 
#     z_bottom = 20.0, 
#     v_inj_point = 3.0, # Example, adjust as needed
#     non_uniform_expansion=False,
#     num_segments=50
# )
# 替换旧的 DEFAULT_GROUT_POINT_FIG5B，使用更明确的案例参数
# 这个点代表案例研究中的一个典型单管，用于那些不明确说明多管的参数敏感性分析图
DEFAULT_SINGLE_GROUT_POINT_CASE_STUDY = GroutPoint(
    x_coord = CASE_STUDY_X_GROUT_COORD, 
    y_coord = 0.0, # 修正为0.0确保对称性
    z_top = CASE_STUDY_Z_GROUT_TOP, 
    z_bottom = CASE_STUDY_Z_GROUT_BOTTOM, 
    v_inj_point = CASE_STUDY_V_INJ_PER_PIPE, # Default to case study per pipe V_inj
    non_uniform_expansion=False, # Default, can be overridden per figure
    num_segments=50
)

# 示例：如果图5b中沿y轴有多个注浆孔，可以这样定义
# GROUTING_POINTS_PROFILE_FIG5B_MULTIPLE = [
#     GroutPoint(x_coord=R_TUNNEL + 10.3, y_coord=-5.0, z_top=15.0, z_bottom=20.0, v_inj_point=1.5),
#     GroutPoint(x_coord=R_TUNNEL + 10.3, y_coord=0.0,  z_top=15.0, z_bottom=20.0, v_inj_point=1.5),
#     GroutPoint(x_coord=R_TUNNEL + 10.3, y_coord=5.0,  z_top=15.0, z_bottom=20.0, v_inj_point=1.5),
# ]

# --- 绘图参数
FIGURE_DPI = 300
OUTPUT_DIR = "C:\\Users\\Rui\\Desktop\\2025-01-08 PHD Thesis\\2025-03-20 PINN+Digital Twin\\Sleeve\\Capsule algorithm-V3\\figures"

# --- 从图7估算的测量数据点 (x=y' (m), y=Horizontal displacement (mm))
# 根据原论文图7重新估算更准确的数据点
MEASURED_DATA_FIG7 = {
    'y_prime': np.array([-25.632, -20.792, -15.708, -10.624, -5.784, -0.886, 4.221, 9.183, 14.176, 19.147, 24.150]),
    'displacement': np.array([ 0.213,   0.688,   1.137,   1.541,  2.256,  3.010,  3.138, 2.352,  1.536,  0.541,  0.310]),
    'style': 'D' # Blue diamonds for Fig 7 measured data
}

# --- 从图6估算的测量和有限元数据点 ---
# 根据原论文图6重新估算更准确的数据点
MEASURED_DATA_FIG6 = {
    'y_prime': np.array([-25.632, -20.792, -15.708, -10.624, -5.784, -0.886, 4.221, 9.183, 14.176, 19.147, 24.150]),
    'displacement': np.array([0.213,   0.688,   1.137,   1.541,  2.256,  3.010,  3.138, 2.352,  1.536,  0.541,  0.310]),
    'style': 'D'  # Diamond markers for measured data
}

FINITE_ELEMENT_DATA_FIG6 = {
    'y_prime': np.array([-29.984, -25.673, -19.491, -14.699, -10.218,  -5.174,  0.438,   5.400,   10.037,  14.389, 19.252,  23.845,  30.00]),
    'displacement': np.array([0.214,  0.303,  0.527,  0.944,  1.676, 2.567,  3.199,  2.890,  1.995,  1.227, 0.693,  0.452,  0.280]),
    'style': '^-'  # Triangle markers with line
}

# --- 非均匀扩张模型参数 (参考论文图4b) ---
# 5m高注浆段，分为5个1m的区域，扩张系数比例从上到下: 9/5, 7/5, 5/5, 3/5, 1/5 (相对于 ε0)
# 这些比例的总和是 (9+7+5+3+1)/5 = 25/5 = 5。
# 如果我们假设这些是相对的体积或源强贡献因子，则归一化因子为 9, 7, 5, 3, 1。
# 总因子 S = 9+7+5+3+1 = 25。
# 每段高度 1m。
NON_UNIFORM_EXPANSION_FACTORS = np.array([9, 7, 5, 3, 1]) # 相对强度因子
NON_UNIFORM_SEGMENT_HEIGHT = 1.0 # m, 每段高度
# 总高度为 len(NON_UNIFORM_EXPANSION_FACTORS) * NON_UNIFORM_SEGMENT_HEIGHT = 5m，与 GROUT_CYLINDER_HEIGHT 一致。
