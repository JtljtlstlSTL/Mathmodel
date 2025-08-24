import numpy as np
from scipy.optimize import fsolve

def compute_radius(θ):
    """ 计算螺线方程 r = 27.2 + (1.7/(2π)) * θ 的极径 r """
    r = 27.2 + (1.7 / (2 * np.pi)) * θ
    return r

def f(x):
    """ 计算函数 f(x) """
    return 0.5 * x * np.sqrt(x**2 + 1) + 0.5 * np.log(np.abs(x + np.sqrt(x**2 + 1)))

def find_r45_point():
    """ 计算 r=4.5m 时的时刻、角度和速度 """
    # 计算 r=4.5m 对应的角度 θ
    # 螺线方程: r = 27.2 + (1.7/(2π)) * θ = 4.5
    # 解方程: θ = (4.5 - 27.2) * (2π) / 1.7
    θ_r45 = (4.5 - 27.2) * (2 * np.pi) / 1.7
    
    # 计算该角度对应的弧长（即时间t）
    # 弧长积分公式: s = (1.7/(2π)) * [f(32π) - f(θ + 32π)]
    term1 = f(32 * np.pi)
    term2 = f(θ_r45 + 32 * np.pi)
    s = (1.7 / (2 * np.pi)) * (term1 - term2)  # 路径长度
    
    # 由于龙头速度恒为1m/s，时间t等于路径长度s
    t_r45 = s
    
    # 计算速度方向（切线方向）
    # 参数方程求导：dr/dθ = 1.7/(2π)
    dr_dθ = 1.7 / (2 * np.pi)
    
    # 切向量分量
    dx_dθ = dr_dθ * np.cos(θ_r45) - compute_radius(θ_r45) * np.sin(θ_r45)
    dy_dθ = dr_dθ * np.sin(θ_r45) + compute_radius(θ_r45) * np.cos(θ_r45)
    
    # 切向量模长
    magnitude = np.sqrt(dx_dθ**2 + dy_dθ**2)
    
    # 单位切向量（速度方向）
    # 由于龙向内运动，θ减小，所以速度方向与θ增加方向相反
    vx = -dx_dθ / magnitude
    vy = -dy_dθ / magnitude
    
    # 速度大小恒为1m/s
    speed = 1.0
    
    # 计算角度（度）
    angle_deg = np.degrees(θ_r45)
    
    return θ_r45, t_r45, vx, vy, speed, angle_deg

# 主程序
if __name__ == "__main__":
    # 计算 r=4.5m 时的参数
    θ, t, vx, vy, speed, angle_deg = find_r45_point()
    
    # 打印结果
    print("\n龙头把手在 r=4.5m 时的参数:")
    print(f"角度 θ = {θ:.6f} 弧度 (约 {angle_deg:.2f} 度)")
    print(f"时间 t = {t:.6f} 秒")
    print(f"速度向量: vx = {vx:.6f} m/s, vy = {vy:.6f} m/s")
    print(f"速度大小: {speed:.6f} m/s (恒定为1m/s)")
    print(f"速度方向: ({vx:.6f}, {vy:.6f})")
    
    # 验证计算结果
    r_calculated = compute_radius(θ)
    print(f"\n验证: 计算得到的 r = {r_calculated:.6f} m (应为 4.5m)")
    
    # 计算运动到该点所需的总圈数
    total_circles = abs(θ) / (2 * np.pi)
    print(f"运动到该点所需的总圈数: {total_circles:.2f} 圈")
    
    # 计算起点和终点角度差
    start_θ = 0
    end_θ = θ
    angle_diff = abs(start_θ - end_θ)
    print(f"起点到终点的角度差: {angle_diff:.2f} 弧度 ({np.degrees(angle_diff):.2f} 度)")