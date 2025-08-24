import numpy as np
import pandas as pd
from scipy.optimize import fsolve


def compute_unit_vector(vx, vy):
    """计算单位向量"""
    magnitude = np.sqrt(vx**2 + vy**2)
    return vx / magnitude, vy / magnitude


def compute_radius(theta, pitch):
    """计算螺线半径，使用变化的螺距参数"""
    return 16*pitch + (pitch / (2 * np.pi)) * theta


def compute_coordinates(theta, pitch):
    """计算坐标"""
    r = compute_radius(theta, pitch)
    return r * np.cos(theta), r * np.sin(theta)


def compute_hand_position(theta, pitch):
    """计算把手位置"""
    return compute_coordinates(theta, pitch)


def f(x):
    """辅助函数"""
    return 0.5 * x * np.sqrt(x**2 + 1) + 0.5 * np.log(np.abs(x + np.sqrt(x**2 + 1)))


def g(x, t, pitch):
    """时间与角度关系函数，使用变化的螺距参数"""
    term1 = f(x + 32 * np.pi)
    term2 = f(32 * np.pi)
    return (pitch / (2 * np.pi)) * (term2 - term1) - t


def solve_theta_of_dragonhead(t_value, pitch):
    """求解龙头角度，使用变化的螺距参数"""
    initial_guess = np.radians(0)
    x_value = fsolve(g, x0=initial_guess, args=(t_value, pitch))
    return x_value[0]


def g_function(theta2, theta1, b, pitch):
    """距离关系函数，使用变化的螺距参数"""
    r1 = compute_radius(theta1, pitch)
    r2 = compute_radius(theta2, pitch)
    return r2**2 - 2 * r1 * r2 * np.cos(theta2 - theta1) + r1**2 - b**2


def compute_theta_next(theta1, b, pitch):
    """计算下一个角度，使用变化的螺距参数"""
    initial_guesses = [theta1 + delta for delta in [0.2, 0.4, 0.7, 0.9]]
    solutions = []
    for guess in initial_guesses:
        theta2_solution = fsolve(g_function, guess, args=(theta1, b, pitch))
        theta2 = theta2_solution[0]
        if 0 < (theta2 - theta1) < 2 * np.pi:
            solutions.append(theta2)
    if solutions:
        return min(np.unique(solutions))
    else:
        return None


def compute_radius_derivative(theta, pitch):
    """计算半径导数，使用变化的螺距参数"""
    return (pitch / (2 * np.pi))


def compute_velocity(theta, pitch):
    """计算速度，使用变化的螺距参数"""
    r = compute_radius(theta, pitch)
    dr = compute_radius_derivative(theta, pitch)
    vx = - r * np.sin(theta) + dr * np.cos(theta)
    vy = r * np.cos(theta) + dr * np.sin(theta)
    return vx, vy


def compute_bench_direction_vector(theta_n, theta_n1, pitch):
    """计算板凳方向向量，使用变化的螺距参数"""
    r_n = compute_radius(theta_n, pitch)
    r_n1 = compute_radius(theta_n1, pitch)
    lx = r_n * np.cos(theta_n) - r_n1 * np.cos(theta_n1)
    ly = r_n * np.sin(theta_n) - r_n1 * np.sin(theta_n1)
    return lx, ly


def compute_velocity_next(theta_n, theta_n1, u_n, pitch):
    """计算下一个速度，使用变化的螺距参数"""
    v_n = compute_unit_vector(*compute_velocity(theta_n, pitch))
    v_n1 = compute_unit_vector(*compute_velocity(theta_n1, pitch))
    l_n_x, l_n_y = compute_bench_direction_vector(theta_n, theta_n1, pitch)
    l_n = np.array([l_n_x, l_n_y])
    numerator = np.dot(v_n, l_n)
    denominator = np.dot(v_n1, l_n)
    if denominator == 0:
        raise ValueError("Denominator is zero, cannot divide.")
    return (numerator / denominator) * u_n


def compute_perpendicular_unit_vector(l_n):
    """计算垂直单位向量"""
    l_n = np.array(l_n)
    h_n = np.array([-l_n[1], l_n[0]])
    return h_n / np.linalg.norm(h_n)


def compute_head_vertices(theta, theta1, pitch):
    """计算龙头顶点，使用变化的螺距参数"""
    r = compute_radius(theta, pitch)
    x0 = r * np.cos(theta)
    y0 = r * np.sin(theta)
    l_0_x, l_0_y = compute_bench_direction_vector(theta, theta1, pitch)
    l0 = np.array(compute_unit_vector(l_0_x, l_0_y))
    h0 = compute_perpendicular_unit_vector(l0)
    center = np.array([x0, y0])
    return [
        center + 0.15 * h0 + 0.275 * l0,
        center + 0.15 * h0 - 3.135 * l0,
        center - 0.15 * h0 + 0.275 * l0,
        center - 0.15 * h0 - 3.135 * l0
    ], l0, h0


def compute_body_vertices(theta, theta1, pitch):
    """计算龙身顶点，使用变化的螺距参数"""
    r = compute_radius(theta, pitch)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    l_n_x, l_n_y = compute_bench_direction_vector(theta, theta1, pitch)
    l_n = np.array(compute_unit_vector(l_n_x, l_n_y))
    h_n = compute_perpendicular_unit_vector(l_n)
    center = np.array([x, y])
    return [
        center + 0.15 * h_n + 0.275 * l_n,
        center + 0.15 * h_n - 1.925 * l_n,
        center - 0.15 * h_n + 0.275 * l_n,
        center - 0.15 * h_n - 1.925 * l_n
    ], l_n, h_n


def projection_interval(vertices, axis):
    """计算投影区间"""
    projections = [np.dot(vertex, axis) for vertex in vertices]
    return min(projections), max(projections)


def rectangles_intersect(vertices1, vertices2):
    """判断两个矩形是否相交"""
    edges1 = [vertices1[i] - vertices1[(i + 1) % 4] for i in range(2)]
    edges2 = [vertices2[i] - vertices2[(i + 1) % 4] for i in range(2)]
    axes = [compute_perpendicular_unit_vector(
        edge) for edge in edges1 + edges2]
    for axis in axes:
        min1, max1 = projection_interval(vertices1, axis)
        min2, max2 = projection_interval(vertices2, axis)
        if max1 < min2 or max2 < min1:
            return False
    return True


def find_min_radius_reached(pitch, max_time=1000):
    """找到龙头能够到达的最小半径（发生碰撞时的半径，或者能盘入的最小半径）"""
    time_array = np.arange(max_time)

    # 计算龙头角度序列
    theta_head_values = np.array(
        [solve_theta_of_dragonhead(t, pitch) for t in time_array])

    # 计算龙身各部分角度
    length_of_dragon_head = 2.86
    distance_of_dragon_body = 1.65
    theta_body_values = np.empty((max_time, 13))

    for t in range(max_time):
        temp = theta_head_values[t]
        for i in range(13):
            if i == 0:
                theta_body_values[t, i] = compute_theta_next(
                    temp, length_of_dragon_head, pitch)
            else:
                theta_body_values[t, i] = compute_theta_next(
                    temp, distance_of_dragon_body, pitch)
            temp = theta_body_values[t, i]

    theta_values_of_each_time = np.hstack(
        (theta_head_values.reshape(-1, 1), theta_body_values))

    # 碰撞检测
    for t in range(max_time):
        theta_0 = theta_values_of_each_time[t, 0]
        theta_1 = theta_values_of_each_time[t, 1]
        vertices_0, l_0, h_0 = compute_head_vertices(theta_0, theta_1, pitch)

        for i in range(2, 13):
            theta_n = theta_values_of_each_time[t, i]
            theta_n1 = theta_values_of_each_time[t, i+1]
            vertices_n, l_n, h_n = compute_body_vertices(
                theta_n, theta_n1, pitch)

            if rectangles_intersect(vertices_0, vertices_n):
                # 发生碰撞，返回此时龙头的半径（这是能到达的最小半径）
                collision_radius = compute_radius(theta_0, pitch)
                return collision_radius

    # 没有发生碰撞，龙头理论上可以盘入到螺线中心
    # 但实际上受时间限制，返回最后时刻的半径
    final_radius = compute_radius(theta_values_of_each_time[-1, 0], pitch)
    return final_radius


def find_pitch_for_target_radius(target_radius=4.5, tolerance=0.001, max_iterations=50):
    """直接求解当龙头到达目标半径时发生碰撞的螺距"""

    def pitch_objective(pitch):
        """目标函数：当使用此螺距时，碰撞半径与目标半径的差值"""
        try:
            collision_radius = find_min_radius_reached(pitch, max_time=2000)
            return collision_radius - target_radius
        except:
            return float('inf')

    # 使用二分法求解
    pitch_min = 0.3
    pitch_max = 0.55

    print(f"直接求解方法：寻找使得碰撞半径等于 {target_radius} 的螺距")
    print("正在计算...")

    for iteration in range(max_iterations):
        pitch_mid = (pitch_min + pitch_max) / 2

        diff = pitch_objective(pitch_mid)

        if abs(diff) < 0.01:  # 精度要求：碰撞半径与目标半径差值小于0.01m
            return pitch_mid

        if diff < 0:
            # 碰撞半径小于目标半径，需要减小螺距
            pitch_max = pitch_mid
        else:
            # 碰撞半径大于目标半径，需要增大螺距
            pitch_min = pitch_mid

    return pitch_max


def binary_search_min_pitch(target_radius=4.5, tolerance=0.001, max_iterations=50):
    """使用二分法寻找最小螺距"""
    # 螺距搜索范围
    pitch_min = 0.3  # 接近0但不能为0
    pitch_max = 0.55   # 原问题二中的螺距

    print(f"二分法搜索，目标半径: {target_radius}")
    print("正在计算...")

    for iteration in range(max_iterations):
        pitch_mid = (pitch_min + pitch_max) / 2

        try:
            min_radius = find_min_radius_reached(pitch_mid)

            if min_radius < target_radius:
                # 能到达的最小半径小于目标半径，说明螺距太小，龙盘入过深，需要更大的螺距
                pitch_max = pitch_mid
            else:
                # 能到达的最小半径大于等于目标半径，说明螺距足够大，可以尝试更小的螺距
                pitch_min = pitch_mid

            # 检查收敛条件
            if abs(pitch_max - pitch_min) < tolerance:
                return pitch_max

        except Exception as e:
            # 出错时尝试减小螺距
            pitch_max = pitch_mid

    print(f"达到最大迭代次数，最小螺距: {pitch_max:.6f}")
    return pitch_max


def main():
    """主函数"""
    print("问题三：寻找最小螺距，使得龙头前把手能够沿着相应的螺线盘入到调头空间的边界")
    print("调头空间：以螺线中心为圆心、直径为9m的圆形区域（半径4.5m）")
    print("=" * 80)

    # 方法1：直接求解当龙头到达目标半径时发生碰撞的螺距
    print("\n方法1：直接反推螺距")
    min_pitch_direct = find_pitch_for_target_radius(target_radius=4.5)

    print("\n" + "=" * 80)
    print(f"方法1结果：最小螺距 = {min_pitch_direct:.6f} m")

    # 验证结果
    print("\n验证方法1结果：")
    final_collision_radius = find_min_radius_reached(min_pitch_direct)
    print(f"使用螺距 {min_pitch_direct:.6f} 时，碰撞半径 = {final_collision_radius:.4f}")

    # 方法2：二分法搜索（原方法，作为对比）
    print("\n" + "=" * 80)
    print("\n方法2：二分法搜索（原方法，作为对比）")
    min_pitch_binary = binary_search_min_pitch(target_radius=4.5)

    print("\n" + "=" * 80)
    print(f"方法2结果：最小螺距 = {min_pitch_binary:.6f} m")

    # 验证结果
    print("\n验证方法2结果：")
    final_collision_radius2 = find_min_radius_reached(min_pitch_binary)
    print(f"使用螺距 {min_pitch_binary:.6f} 时，碰撞半径 = {final_collision_radius2:.4f}")

    print("\n" + "=" * 80)
    print("最终对比：")
    print(
        f"方法1（直接反推）：螺距 = {min_pitch_direct:.6f} m，碰撞半径 = {final_collision_radius:.4f} m")
    print(
        f"方法2（二分搜索）：螺距 = {min_pitch_binary:.6f} m，碰撞半径 = {final_collision_radius2:.4f} m")

    # 选择更准确的结果
    if abs(final_collision_radius - 4.5) <= abs(final_collision_radius2 - 4.5):
        best_pitch = min_pitch_direct
        best_radius = final_collision_radius
        best_method = "直接反推"
    else:
        best_pitch = min_pitch_binary
        best_radius = final_collision_radius2
        best_method = "二分搜索"

    print(f"\n推荐结果（{best_method}）：螺距 = {best_pitch:.6f} m")

    if best_radius <= 4.5:
        print("✓ 验证成功：龙头前把手能够到达调头空间边界")
    else:
        print("✗ 验证失败：需要调整螺距")

    # 保存结果
    with open('result_q3_min_pitch.txt', 'w', encoding='utf-8') as f:
        f.write(f"问题三结果：\n")
        f.write(f"最小螺距（直接反推）: {min_pitch_direct:.6f} m\n")
        f.write(f"碰撞半径（直接反推）: {final_collision_radius:.4f} m\n")
        f.write(f"最小螺距（二分搜索）: {min_pitch_binary:.6f} m\n")
        f.write(f"碰撞半径（二分搜索）: {final_collision_radius2:.4f} m\n")
        f.write(f"推荐结果: {best_pitch:.6f} m ({best_method})\n")
        f.write(f"调头空间半径: 4.5 m\n")
        f.write(f"是否满足条件: {'是' if best_radius <= 4.5 else '否'}\n")

    print(f"\n结果已保存到 'result_q3_min_pitch.txt'")


if __name__ == "__main__":
    main()
