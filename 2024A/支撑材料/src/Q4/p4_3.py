import numpy as np
from scipy.optimize import fsolve

def compute_radius(θ):
    """ 计算螺线方程 r = 27.2 + (1.7/(2π)) * θ 的极径 r """
    r = 27.2 + (1.7 / (2 * np.pi)) * θ
    return r

def compute_coordinates(θ):
    """ 计算点的坐标 (r cos(θ), r sin(θ)) """
    r = compute_radius(θ)
    x = r * np.cos(θ)
    y = r * np.sin(θ)
    return (x, y)

def f(x):
    """ 计算函数 f(x) """
    return 0.5 * x * np.sqrt(x**2 + 1) + 0.5 * np.log(np.abs(x + np.sqrt(x**2 + 1)))

def g(x, t):
    """ 计算函数 g(x) - t """
    term1 = f(x + 32 * np.pi)
    term2 = f(32 * np.pi)
    return (1.7 / (2 * np.pi)) * (term2 - term1) - t

def solve_θ_of_dragonhead(t_value):
    """ 解方程 g(x) = t 的 x 值 """
    initial_guess = np.radians(0)
    try:
        x_value = fsolve(g, x0=initial_guess, args=(t_value))
        return x_value[0]
    except Exception as e:
        print(f"Failed to solve: {e}")
        return None

def compute_velocity_components(θ):
    """ 计算速度分量 """
    r = compute_radius(θ)
    dr_dθ = 1.7 / (2 * np.pi)  # 螺距为1.7m
    
    # 切向量分量
    dx_dθ = dr_dθ * np.cos(θ) - r * np.sin(θ)
    dy_dθ = dr_dθ * np.sin(θ) + r * np.cos(θ)
    
    # 切向量模长
    magnitude = np.sqrt(dx_dθ**2 + dy_dθ**2)
    
    # 单位切向量（速度方向）
    # 由于龙向内运动，θ减小，所以速度方向与θ增加方向相反
    vx = -dx_dθ / magnitude
    vy = -dy_dθ / magnitude
    
    # 速度大小恒为1m/s
    speed = 1.0
    
    return vx, vy, speed

def compute_velocity_components_with_speed(θ, prev_tangent=None, prev_speed=1.0):
    """ 计算速度分量和大小，支持链式递推 """
    r = compute_radius(θ)
    dr_dθ = 1.7 / (2 * np.pi)
    dx_dθ = dr_dθ * np.cos(θ) - r * np.sin(θ)
    dy_dθ = dr_dθ * np.sin(θ) + r * np.cos(θ)
    magnitude = np.sqrt(dx_dθ**2 + dy_dθ**2)
    tangent = np.array([-dx_dθ / magnitude, -dy_dθ / magnitude])  # 速度方向
    if prev_tangent is not None:
        # 链式速度递推
        cos_angle = np.clip(np.dot(tangent, prev_tangent), -1.0, 1.0)
        speed = prev_speed * cos_angle
        speed = np.abs(speed)  # 速度大小不能为负
    else:
        speed = 1.0  # 龙头速度
    vx = tangent[0] * speed
    vy = tangent[1] * speed
    return vx, vy, speed, tangent


def compute_body_point(θ_head, section_index, body_length=1.65):
    """ 计算龙身点的位置和速度 """
    # 计算龙身点的角度（简化模型：沿切线方向延伸）
    # 实际应用中可能需要更复杂的模型
    θ_body = θ_head - (section_index + 1) * 0.01  # 简化处理
    
    # 计算位置
    x, y = compute_coordinates(θ_body)
    
    # 计算速度
    vx, vy, speed = compute_velocity_components(θ_body)
    
    return (x, y, vx, vy, speed)

def compute_body_points_chain(θ_head, sections):
    """链式递推计算所有龙身点的位置和速度"""
    points = {}
    θ_prev = θ_head
    vx_prev, vy_prev, speed_prev, tangent_prev = compute_velocity_components_with_speed(θ_head)
    for idx, section_index in enumerate(sections):
        θ_body = θ_head - (section_index + 1) * 0.01
        vx, vy, speed, tangent = compute_velocity_components_with_speed(θ_body, tangent_prev, speed_prev)
        x, y = compute_coordinates(θ_body)
        points[section_index] = (x, y, vx, vy, speed)
        # 为下一个点递推
        θ_prev = θ_body
        vx_prev, vy_prev, speed_prev, tangent_prev = vx, vy, speed, tangent
    return points

def compute_tail_point(θ_head):
    """ 计算龙尾后把手的位置和速度 """
    # 计算龙尾的角度（简化模型）
    θ_tail = θ_head - 223 * 0.01  # 223个龙身点
    
    # 计算位置
    x, y = compute_coordinates(θ_tail)
    
    # 计算速度
    vx, vy, speed = compute_velocity_components(θ_tail)
    
    return (x, y, vx, vy, speed)

def compute_tail_point_chain(θ_head, last_tangent, last_speed):
    θ_tail = θ_head - 223 * 0.01
    vx, vy, speed, tangent = compute_velocity_components_with_speed(θ_tail, last_tangent, last_speed)
    x, y = compute_coordinates(θ_tail)
    return (x, y, vx, vy, speed)


def compute_dragon_state(t):
    """ 计算给定时间t的龙各部分状态（链式速度递推） """
    # 计算龙头角度
    θ_head = solve_θ_of_dragonhead(t)
    
    # 计算龙头位置
    x_head, y_head = compute_coordinates(θ_head)
    
    # 计算龙头速度
    vx_head, vy_head, speed_head, tangent_head = compute_velocity_components_with_speed(θ_head)
    
    # 计算特定龙身点的位置和速度
    sections = [1, 51, 101, 151, 201]
    body_points = {}
    θ_prev = θ_head
    vx_prev, vy_prev, speed_prev, tangent_prev = vx_head, vy_head, speed_head, tangent_head
    for section_index in sections:
        θ_body = θ_head - (section_index + 1) * 0.01
        vx, vy, speed, tangent = compute_velocity_components_with_speed(θ_body, tangent_prev, speed_prev)
        x, y = compute_coordinates(θ_body)
        body_points[section_index] = (x, y, vx, vy, speed)
        θ_prev = θ_body
        vx_prev, vy_prev, speed_prev, tangent_prev = vx, vy, speed, tangent
    # 龙尾
    x_tail, y_tail, vx_tail, vy_tail, speed_tail = compute_tail_point_chain(θ_head, tangent_prev, speed_prev)
    return {
        "time": t,
        "head": {
            "position": (x_head, y_head),
            "velocity": (vx_head, vy_head, speed_head)
        },
        "body": body_points,
        "tail": {
            "position": (x_tail, y_tail),
            "velocity": (vx_tail, vy_tail, speed_tail)
        }
    }

def format_result(result):
    """ 格式化结果输出 """
    t = result["time"]
    output = f"时间: {t:.6f} 秒\n"
    
    # 龙头前把手
    head_pos = result["head"]["position"]
    head_vel = result["head"]["velocity"]
    output += "龙头前把手:\n"
    output += f"  位置: ({head_pos[0]:.6f}, {head_pos[1]:.6f}) m\n"
    output += f"  速度: ({head_vel[0]:.6f}, {head_vel[1]:.6f}) m/s, 大小: {head_vel[2]:.6f} m/s\n"
    
    # 特定龙身点
    for section in [1, 51, 101, 151, 201]:
        data = result["body"][section]
        output += f"龙头后面第{section}节龙身前把手:\n"
        output += f"  位置: ({data[0]:.6f}, {data[1]:.6f}) m\n"
        output += f"  速度: ({data[2]:.6f}, {data[3]:.6f}) m/s, 大小: {data[4]:.6f} m/s\n"
    
    # 龙尾后把手
    tail_pos = result["tail"]["position"]
    tail_vel = result["tail"]["velocity"]
    output += "龙尾后把手:\n"
    output += f"  位置: ({tail_pos[0]:.6f}, {tail_pos[1]:.6f}) m\n"
    output += f"  速度: ({tail_vel[0]:.6f}, {tail_vel[1]:.6f}) m/s, 大小: {tail_vel[2]:.6f} m/s\n"
    
    return output

# 主程序
if __name__ == "__main__":
    # 定义时间点
    times = [1230.042542, 1280.042542, 1330.042542]
    
    # 计算并输出结果
    for t in times:
        result = compute_dragon_state(t)
        print(format_result(result))
        print("-" * 80)