D = 120; % 海水深度（米）

% 测线位置（单位：海里）
y_nm = [0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1];
y = y_nm * 1852; % 转换为米

% 初始坡度角和开角
alpha0_deg = 1.5;
theta_deg = 120;

% 转为弧度
alpha0 = deg2rad(alpha0_deg);
theta = deg2rad(theta_deg);

% 方向角 beta（单位：度）
beta_deg = [0 45 90 135 180 225 270 315];

% 初始化二维矩阵（行：y，列：β）
W = zeros(length(y), length(beta_deg)); % 覆盖宽度
H = zeros(length(y), length(beta_deg)); % 水深

% 二维循环计算 W(i,j) = w(yᵢ, βⱼ)
for i = 1 : length(y)
    for j = 1 : length(beta_deg)
        beta = deg2rad(beta_deg(j));
        x_rad = atan(-cos(beta) * tan(alpha0));
        h = D - y(i) * tan(x_rad);
        x_rad = atan2(tan(alpha0) * sin(beta), 1);
        w = h * (1 / cos(x_rad + theta / 2) + 1 / cos(x_rad - theta / 2)) * sin(theta / 2);
        H(i, j) = h;
        W(i, j) = w;
    end
end

W_swapped = W';
H_swapped = H';

% 打印标题
fprintf('=== 覆盖宽度 w(beta, y) 表（单位：米，保留4位小数） ===\n');

% 打印列标题（y）
fprintf('%10s', '');  % 左上角空位
for i = 1:length(y_nm)
    fprintf('%12s', sprintf('y=%.1fNM', y_nm(i)));
end
fprintf('\n');

% 打印每一行（beta 对应的行）
for j = 1:length(beta_deg)
    fprintf('%10s', sprintf('Beta_%d', beta_deg(j)));
    for i = 1:length(y_nm)
        fprintf('%12.4f', W_swapped(j, i));  % 保留4位小数
    end
    fprintf('\n');
end
