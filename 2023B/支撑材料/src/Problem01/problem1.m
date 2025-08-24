% 输入参数
D = 70;  % 海水深度
y = [-800 -600 -400 -200 0 200 400 600 800];  % 测线位置
alpha_deg = 1.5;  % 坡度角，单位：度
theta_deg = 120;  % 开角，单位：度

% 调用函数计算
[h, w, t, t_str] = calculate_overlap_rate(D, y, alpha_deg, theta_deg);