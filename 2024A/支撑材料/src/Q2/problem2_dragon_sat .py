import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# === 参数设置 ===
def compute_unit_vector(vx, vy):
    magnitude = np.sqrt(vx**2 + vy**2)
    return vx / magnitude, vy / magnitude

def compute_radius(theta):
    return 8.8 + (0.55 / (2 * np.pi)) * theta

def compute_coordinates(theta):
    r = compute_radius(theta)
    return r * np.cos(theta), r * np.sin(theta)

def compute_hand_position(theta):
    return compute_coordinates(theta)

def f(x):
    return 0.5 * x * np.sqrt(x**2 + 1) + 0.5 * np.log(np.abs(x + np.sqrt(x**2 + 1)))

def g(x, t):
    term1 = f(x + 32 * np.pi)
    term2 = f(32 * np.pi)
    return (0.55 / (2 * np.pi)) * (term2 - term1) - t

def solve_theta_of_dragonhead(t_value):
    initial_guess = np.radians(0)
    x_value = fsolve(g, x0=initial_guess, args=(t_value))
    return x_value[0]

def g_function(theta2, theta1, b):
    r1 = compute_radius(theta1)
    r2 = compute_radius(theta2)
    return r2**2 - 2 * r1 * r2 * np.cos(theta2 - theta1) + r1**2 - b**2

def compute_theta_next(theta1, b):
    initial_guesses = [theta1 + delta for delta in [0.2, 0.4, 0.7, 0.9]]
    solutions = []
    for guess in initial_guesses:
        theta2_solution = fsolve(g_function, guess, args=(theta1, b))
        theta2 = theta2_solution[0]
        if 0 < (theta2 - theta1) < 2 * np.pi:
            solutions.append(theta2)
    if solutions:
        return min(np.unique(solutions))
    else:
        return None

def compute_radius_derivative(theta):
    return (0.55 / (2 * np.pi))

def compute_velocity(theta):
    r = compute_radius(theta)
    dr = compute_radius_derivative(theta)
    vx = - r * np.sin(theta) + dr * np.cos(theta)
    vy = r * np.cos(theta) + dr * np.sin(theta)
    return vx, vy

def compute_bench_direction_vector(theta_n, theta_n1):
    r_n = compute_radius(theta_n)
    r_n1 = compute_radius(theta_n1)
    lx = r_n * np.cos(theta_n) - r_n1 * np.cos(theta_n1)
    ly = r_n * np.sin(theta_n) - r_n1 * np.sin(theta_n1)
    return lx, ly

def compute_velocity_next(theta_n, theta_n1, u_n):
    v_n = compute_unit_vector(*compute_velocity(theta_n))
    v_n1 = compute_unit_vector(*compute_velocity(theta_n1))
    l_n_x, l_n_y = compute_bench_direction_vector(theta_n, theta_n1)
    l_n = np.array([l_n_x, l_n_y])
    numerator = np.dot(v_n, l_n)
    denominator = np.dot(v_n1, l_n)
    if denominator == 0:
        raise ValueError("Denominator is zero, cannot divide.")
    return (numerator / denominator) * u_n

def compute_perpendicular_unit_vector(l_n):
    l_n = np.array(l_n)
    h_n = np.array([-l_n[1], l_n[0]])
    return h_n / np.linalg.norm(h_n)

def compute_head_vertices(theta, theta1):
    r = compute_radius(theta)
    x0 = r * np.cos(theta)
    y0 = r * np.sin(theta)
    l_0_x, l_0_y = compute_bench_direction_vector(theta, theta1)
    l0 = np.array(compute_unit_vector(l_0_x, l_0_y))
    h0 = compute_perpendicular_unit_vector(l0)
    center = np.array([x0, y0])
    return [
        center + 0.15 * h0 + 0.275 * l0,
        center + 0.15 * h0 - 3.135 * l0,
        center - 0.15 * h0 + 0.275 * l0,
        center - 0.15 * h0 - 3.135 * l0
    ], l0, h0

def compute_body_vertices(theta, theta1):
    r = compute_radius(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    l_n_x, l_n_y = compute_bench_direction_vector(theta, theta1)
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
    projections = [np.dot(vertex, axis) for vertex in vertices]
    return min(projections), max(projections)

def rectangles_intersect(vertices1, vertices2):
    edges1 = [vertices1[i] - vertices1[(i + 1) % 4] for i in range(2)]
    edges2 = [vertices2[i] - vertices2[(i + 1) % 4] for i in range(2)]
    axes = [compute_perpendicular_unit_vector(edge) for edge in edges1 + edges2]
    for axis in axes:
        min1, max1 = projection_interval(vertices1, axis)
        min2, max2 = projection_interval(vertices2, axis)
        if max1 < min2 or max2 < min1:
            return False
    return True

# 检测是否发生碰撞

def detect_collision_at_time(t_float):
    theta_head = solve_theta_of_dragonhead(t_float)
    temp = theta_head
    theta_segments = []
    for i in range(13):
        b = 2.86 if i == 0 else 1.65
        theta_next = compute_theta_next(temp, b)
        theta_segments.append(theta_next)
        temp = theta_next
    thetas = [theta_head] + theta_segments
    theta_0 = thetas[0]
    theta_1 = thetas[1]
    vertices_0, _, _ = compute_head_vertices(theta_0, theta_1)
    for i in range(2, 13):
        theta_n = thetas[i]
        theta_n1 = thetas[i+1]
        vertices_n, _, _ = compute_body_vertices(theta_n, theta_n1)
        if rectangles_intersect(vertices_0, vertices_n):
            return True
    return False

# 二分查找精确碰撞时间

def find_precise_collision_time(t_low=411.0, t_high=413.0, tol=1e-6):
    while t_high - t_low > tol:
        mid = (t_low + t_high) / 2
        if detect_collision_at_time(mid):
            t_high = mid
        else:
            t_low = mid
    return (t_low + t_high) / 2

print("开始精确检测碰撞时间...")
t_precise = find_precise_collision_time()
print(f"更精确的碰撞发生时间为 t = {t_precise:.6f} 秒")
