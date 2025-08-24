import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# 1. 参数设置（根据题目精确调整）
initial_theta = 32 * np.pi  # 初始在第16圈（2π*16=32π）
pitch = 0.55                # 螺距0.55米
head_handle_dist = 2.86     # 龙头把手间距2.86米（题目图3）
body_handle_dist = 1.65     # 龙身把手间距1.65米（题目说明）
total_sections = 223        # 总节数: 1龙头 + 221龙身 + 1龙尾

# 2. 螺线方程修正（关键修正点）
def compute_radius(theta):
    """r = r0 + (h/2π) * (θ - θ0)"""
    return 8.8 + (pitch / (2 * np.pi)) * (theta - initial_theta)

# 3. 坐标计算
def compute_coordinates(theta):
    r = compute_radius(theta)
    x = np.round(r * np.cos(theta), 6)
    y = np.round(r * np.sin(theta), 6)
    return x, y

# 4. 时间-角度关系（完全重写）
def theta_at_time(t):
    """数值求解θ(t)满足 ∫√(r² + (dr/dθ)²)dθ = vt"""
    def equation(theta):
        integral = (pitch/(2*np.pi)) * (
            0.5*(theta*np.sqrt(theta**2 + 1) + np.log(theta + np.sqrt(theta**2 + 1))) -
            0.5*(initial_theta*np.sqrt(initial_theta**2 + 1) + np.log(initial_theta + np.sqrt(initial_theta**2 + 1)))
        )
        return integral - t
    return fsolve(equation, x0=initial_theta - t/10)[0]

# 5. 计算龙头轨迹（0-300秒）
time_points = [0, 60, 120, 180, 240, 300]
theta_head = [theta_at_time(t) for t in time_points]

# 6. 计算龙身和龙尾轨迹（精确几何约束）
theta_body = np.zeros((len(time_points), total_sections-1))

def find_next_theta(theta1, L):
    """精确求解相邻点θ2满足距离约束"""
    r1 = compute_radius(theta1)
    def equation(theta2):
        r2 = compute_radius(theta2)
        return r2**2 + r1**2 - 2*r1*r2*np.cos(theta2-theta1) - L**2
    
    # 使用更稳健的求解方法
    theta_guess = theta1 + L / r1  # 基于弧长的初始猜测
    solution = fsolve(equation, x0=theta_guess, xtol=1e-6)
    return solution[0]

for i, t in enumerate(time_points):
    # 龙头到第一节龙身（特殊距离2.86m）
    theta_body[i, 0] = find_next_theta(theta_head[i], head_handle_dist)
    
    # 后续龙身点（距离1.65m）
    for j in range(1, total_sections-1):
        theta_body[i, j] = find_next_theta(theta_body[i, j-1], body_handle_dist)

# 7. 选择特定点
selected_indices = [0, 1, 51, 101, 151, 201, 221]  # 龙头(0),第1节(1),...,龙尾(221)

# 8. 生成结果表格
result = []
for i, t in enumerate(time_points):
    theta_all = np.concatenate(([theta_head[i]], theta_body[i]))
    theta_selected = theta_all[selected_indices]
    x, y = compute_coordinates(theta_selected)
    result.append((t, x, y))

# 9. 转换为DataFrame
df_result = pd.DataFrame({
    '时间(s)': [r[0] for r in result]*7,
    '点': np.repeat(['龙头', '第1节龙身', '第51节龙身', '第101节龙身', 
                   '第151节龙身', '第201节龙身', '龙尾(后)'], 6),
    'x(m)': np.concatenate([r[1] for r in result]),
    'y(m)': np.concatenate([r[2] for r in result])
})

# 10. 验证输出
print("修正后的龙头初始位置验证：")
print(f"t=0s时龙头坐标: {compute_coordinates(theta_head[0])} (应与表3的(8.8, 0.0)一致)")

# 11. 保存结果
df_result.to_csv('Q1_final_result.csv', index=False)
print("结果已保存到Q1_final_result.csv")