import numpy as np

# === 参数设置 === #
theta_ideal = np.pi * 7 / 18    # 理想角 70°
delta = 0.01                     # 每轮迭代步长
max_iter = 1000                   # 最大迭代轮数

def polar_to_cartesian(r, theta_deg):
    theta_rad = np.radians(theta_deg)
    return r * np.array([np.cos(theta_rad), np.sin(theta_rad)])

def angle_between_vectors(v1, v2):
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.arccos(np.clip(cos_theta, -1, 1))

def normalize_angle(theta):
    return theta % (2 * np.pi)

def get_unit_vector(angle_rad):
    return np.array([np.cos(angle_rad), np.sin(angle_rad)])


#迭代修正函数
def iterate_position_correction(initial_positions, center=np.array([0.0, 0.0]), target_r=None):
    positions = initial_positions.copy()
    n = len(positions)
    if target_r is None:
        # 默认目标半径为初始极径均值
        target_r = np.mean([np.linalg.norm(p - center) for p in positions])

    for iteration in range(max_iter):
        print(f"\n 第 {iteration + 1} 轮迭代")
        new_positions = positions.copy()

        for i in range(n):  # 所有无人机都参与
            Pi = positions[i]
            Pi_prev = positions[(i - 1) % n]
            Pi_next = positions[(i + 1) % n]

            # 向量定义
            vec_prev = Pi_prev - Pi
            vec_center = center - Pi
            vec_next = Pi_next - Pi

            # 两个夹角
            theta1 = angle_between_vectors(vec_prev, vec_center)
            theta2 = angle_between_vectors(vec_next, vec_center)

            # 偏离程度（角度偏差归一化）
            d1 = abs(theta1 - theta_ideal)
            d2 = abs(theta2 - theta_ideal)
            norm = np.sqrt(d1**2 + d2**2) + 1e-6
            c1, c2 = d1 / norm, d2 / norm

            # 两个角平分线方向（修正）
            bisector1 = (vec_prev / np.linalg.norm(vec_prev) + vec_center / np.linalg.norm(vec_center))
            bisector1 /= np.linalg.norm(bisector1)
            bisector2 = (vec_next / np.linalg.norm(vec_next) + vec_center / np.linalg.norm(vec_center))
            bisector2 /= np.linalg.norm(bisector2)
            direction = c1 * bisector1 + c2 * bisector2
            direction /= np.linalg.norm(direction)

            # 径向修正分量
            r_now = np.linalg.norm(Pi - center)
            radial_dir = (Pi - center) / (r_now + 1e-8)
            radial_correction = (target_r - r_now) * 0.2 * radial_dir  # 0.2为径向收敛速率

            move_vec = delta * direction + radial_correction
            new_positions[i] += move_vec

            print(f"FY0{i+1}: θ1={np.degrees(theta1):.2f}°, θ2={np.degrees(theta2):.2f}°, dx={move_vec[0]:.3f}, dy={move_vec[1]:.3f}, r={np.linalg.norm(new_positions[i]-center):.2f}")

        positions = new_positions.copy()

    return positions


# === 主程序入口 === #

if __name__ == "__main__":
    # 初始极坐标位置（单位：米 / 度）
    init_polar = {
        1: (100, 0),
        2: (98, 40.10),
        3: (112, 80.21),
        4: (105, 119.75),
        5: (98, 159.86),
        6: (112, 199.96),
        7: (105, 240.07),
        8: (98, 280.17),
        9: (112, 320.28)
    }

    # 转换为直角坐标
    initial_positions = [polar_to_cartesian(r, theta) for r, theta in init_polar.values()]

    # 执行迭代修正
    corrected_positions = iterate_position_correction(initial_positions)

    # 输出最终位置（极坐标形式）
    print("\n✅ 最终修正位置坐标：")
    for i, pos in enumerate(corrected_positions):
        r = np.linalg.norm(pos)
        theta = np.degrees(np.arctan2(pos[1], pos[0])) % 360
        print(f"FY0{i+1}: x={pos[0]:.2f}, y={pos[1]:.2f}, r={r:.2f}, θ={theta:.2f}°")
