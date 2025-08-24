import numpy as np

def locate_uav_correct(R, t, i, j, alpha1_deg, alpha2_deg):
    if i == j:
        raise ValueError("信号源编号 i 和 j 不能相同")
    if not (1 <= i <= 9 and 1 <= j <= 9):
        raise ValueError("编号必须在 1 到 9 之间")

    alpha1 = np.radians(alpha1_deg)
    alpha2 = np.radians(alpha2_deg)
    theta_i = np.radians(40 * i - 40)
    theta_j = np.radians(40 * j - 40)

    # 定义单位圆方向角
    theta_t = np.radians(40 * t - 40)

    # 顺时针方向的夹角（单位弧度）
    angle_cw1 = (theta_t - theta_i) % (2 * np.pi)
    angle_cw2 = (theta_t - theta_j) % (2 * np.pi)

    # 替代原来的 if 判断
    if angle_cw1 < np.pi:
       theta1 = np.pi/2 + theta_i - alpha1
    else:
       theta1 = -np.pi/2 + theta_i + alpha1

    if angle_cw2 < np.pi:
        theta2 = np.pi/2 + theta_j - alpha2
    else:
        theta2 = -np.pi/2 + theta_j + alpha2

    r1 = R / (2 * np.sin(alpha1))
    r2 = R / (2 * np.sin(alpha2))

    numerator = np.cos(theta2) * np.sin(alpha1) - \
        np.cos(theta1) * np.sin(alpha2)
    denominator = np.sin(theta2) * np.sin(alpha1) - \
        np.sin(theta1) * np.sin(alpha2)
    theta0_base = np.arctan2(-numerator, denominator) % (2 * np.pi)

    # π多解消除
    def total_alpha(theta_test):
        # 极角转为坐标点
        B = np.array([np.cos(theta_test), np.sin(theta_test)])  # FY0t
        A1 = np.array([np.cos(theta_i), np.sin(theta_i)])       # FY0i
        A2 = np.array([np.cos(theta_j), np.sin(theta_j)])       # FY0j
        C = np.array([0.0, 0.0])                                 # FY00

    # ∠A1–B–C = α₁
        BA1 = A1 - B
        BC = C - B
        cos1 = np.dot(BA1, BC) / (np.linalg.norm(BA1) * np.linalg.norm(BC))
        a1 = np.arccos(np.clip(cos1, -1, 1))

    # ∠A2–B–C = α₂
        BA2 = A2 - B
        cos2 = np.dot(BA2, BC) / (np.linalg.norm(BA2) * np.linalg.norm(BC))
        a2 = np.arccos(np.clip(cos2, -1, 1))

        return a1 + a2

    expected_sum = alpha1 + alpha2
    candidates = [theta0_base, (theta0_base + np.pi) % (2 * np.pi)]
    theta0 = min(candidates, key=lambda th: abs(
        total_alpha(th) - expected_sum))

    cos_val = np.cos(theta0 - theta1)
    r0 = 2 * r1 * cos_val

    theta0_deg = np.degrees(theta0)
    return r0, theta0_deg


# === 主程序入口 ===
if __name__ == "__main__":
    R = 100
    t = 2
    i = 1
    j = 3
    alpha1 = 70
    alpha2 = 70

    r0, theta0 = locate_uav_correct(R, t, i, j, alpha1, alpha2)
    print(f"✅ 定位结果：r₀ = {r0:.4f} m, θ₀ = {theta0:.4f}°")
