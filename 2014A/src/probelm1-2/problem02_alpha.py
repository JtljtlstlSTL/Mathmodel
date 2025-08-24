import numpy as np

def compute_angle_FY00_T_I(T, I):
    if T == I:
        raise ValueError("T 和 I 不能相同")

    # 每个点的角度（以 40° 间隔，转换为弧度）
    theta_T = np.radians(40 * T - 40)
    theta_I = np.radians(40 * I - 40)

    # 极坐标单位向量表示
    P0 = np.array([0.0, 0.0])  # 原点 FY00
    PT = np.array([np.cos(theta_T), np.sin(theta_T)])  # FY0T
    PI = np.array([np.cos(theta_I), np.sin(theta_I)])  # FY0I

    # ∠P0 - PT - PI
    vec1 = P0 - PT
    vec2 = PI - PT

    dot_product = np.dot(vec1, vec2)
    norm_product = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    cos_angle = np.clip(dot_product / norm_product, -1.0, 1.0)
    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)

    return angle_deg

if __name__ == "__main__":
    print("角 ∠FY00-FY0T-FY0I（单位：度）：")
    for T in range(2, 10):         # T 从 2 到 9
        for I in range(1, 10):     # I 从 1 到 9
            if I == T:
                continue
            angle = compute_angle_FY00_T_I(T, I)
            print(f"T = {T}, I = {I} → ∠FY00-FY0{T}-FY0{I} = {angle:.4f}°")

    print("角 ∠FY0J-FY0T-FY01（单位：度）：")
    for J in range(2, 10):         # J 从 2 到 9
        for T in range(2, 10):     # T 从 2 到 9
            if J == T:
                continue
            angle = compute_angle_FY00_T_I(J, T)
            print(f"J = {J}, T = {T} → ∠FY0{J}-FY0{T}-FY01 = {angle:.4f}°")
