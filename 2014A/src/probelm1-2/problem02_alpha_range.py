import numpy as np

def compute_angle_range_FY00_T_I(T_idx, I_idx, radius=15, num_samples=360):
    if T_idx == I_idx:
        raise ValueError("T 和 I 不能相同")
    R=100

    # 原始 FY0T 角度位置（单位圆上的方向）
    theta_T = np.radians(40 * T_idx - 40)
    center_T = np.array([R*np.cos(theta_T), R*np.sin(theta_T)])  # 圆心坐标

    # FY0I 位置
    theta_I = np.radians(40 * I_idx - 40)
    point_I = np.array([R*np.cos(theta_I), R*np.sin(theta_I)])  # FY0I
    point_0 = np.array([0.0, 0.0])  # FY00

    # 采样点
    angles = np.linspace(0, 2 * np.pi, num_samples, endpoint=False)
    angle_list = []

    for phi in angles:
        # 当前采样点 T'
        T_sample = center_T + radius * np.array([np.cos(phi), np.sin(phi)])

        # ∠FY00 - T_sample - FY0I
        vec1 = point_0 - T_sample
        vec2 = point_I - T_sample
        dot_product = np.dot(vec1, vec2)
        norm_product = np.linalg.norm(vec1) * np.linalg.norm(vec2)
        cos_angle = np.clip(dot_product / norm_product, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)
        angle_list.append(angle_deg)

    return min(angle_list), max(angle_list)


def compute_angle_range_FY0J_T_1(J_idx, T_idx, radius=15, num_samples=360):
    if T_idx == J_idx:
        raise ValueError("T 和 J 不能相同")
    R=100

    # 原始 FY0T 角度位置（单位圆上的方向）
    theta_T = np.radians(40 * T_idx - 40)
    center_T = np.array([R*np.cos(theta_T), R*np.sin(theta_T)])  # 圆心坐标

    # FY0J 位置
    theta_J = np.radians(40 * J_idx - 40)
    point_J = np.array([R*np.cos(theta_J), R*np.sin(theta_J)])  # FY0J
    point_0 = np.array([0.0, 0.0])  # FY00

    # 采样点
    angles = np.linspace(0, 2 * np.pi, num_samples, endpoint=False)
    angle_list = []

    for phi in angles:
        # 当前采样点 T'
        T_sample = center_T + radius * np.array([np.cos(phi), np.sin(phi)])

        # ∠FY00 - T_sample - FY0J
        vec1 = point_0 - T_sample
        vec2 = point_J - T_sample
        dot_product = np.dot(vec1, vec2)
        norm_product = np.linalg.norm(vec1) * np.linalg.norm(vec2)
        cos_angle = np.clip(dot_product / norm_product, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        angle_deg = np.degrees(angle_rad)
        angle_list.append(angle_deg)

    return min(angle_list), max(angle_list)


if __name__ == "__main__":
    print("角 ∠FY00-FY0T-FY0I（单位：度）：")
    for T in range(2, 10):         # T 从 2 到 9
        for I in range(1, 10):     # I 从 1 到 9
            if I == T:
                continue
            angle_min, angle_max = compute_angle_range_FY00_T_I(T, I)
            print(f"T = {T}, I = {I} → ∠FY00-FY0{T}-FY0{I} = [{angle_min:.4f}°, {angle_max:.4f}°]")

    print("角 ∠FY0J-FY0T-FY01（单位：度）：")
    for J in range(2, 10):         # J 从 2 到 9
        for T in range(2, 10):     # T 从 2 到 9
            if J == T:
                continue
            angle_min, angle_max = compute_angle_range_FY0J_T_1(J, T)
            print(f"J = {J}, T = {T} → ∠FY0{J}-FY0{T}-FY01 = [{angle_min:.4f}°, {angle_max:.4f}°]")
