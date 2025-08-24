clc; clear;

data = xlsread('附件.xlsx');

% 提取坐标向量
x = data(1, 2:end);              % 横坐标（单位：海里）
y = data(2:end, 1);              % 纵坐标（单位：海里）
z_nm = data(2:end, 2:end);       % 深度（单位：米）
z = z_nm / 1852;                 % 转换为海里

% 原始网格
[X, Y] = meshgrid(x, y);

[xi, yi] = meshgrid(linspace(min(x), max(x), 200), ...
                    linspace(min(y), max(y), 200));
zi = griddata(X, Y, z, xi, yi, 'cubic');

figure;
contour(xi, yi, zi, 50, 'k'); hold on;
xlabel('Longitude (NM)');
ylabel('Latitude (NM)');
title('批量平面拟合与等深线叠加');
axis equal;
view(3);
colormap(parula);


regions = [
    0, 1.1, 0, 2;
    1.1, 2.3, 0, 2;
    0, 1.1, 2, 3.9;
    1.1, 2.3, 2, 3.9;
    0, 2.3, 3.9, 5;
    2.3, 4, 0, 3;
    2.3, 4, 3, 4.4;
    2.3, 4, 4.4, 5;
];


for k = 1:size(regions, 1)
    x_range = regions(k, 1:2);
    y_range = regions(k, 3:4);

    % 提取该区域内的点
    X_vec = xi(:); Y_vec = yi(:); Z_vec = zi(:);
    in_region = X_vec >= x_range(1) & X_vec <= x_range(2) & ...
                Y_vec >= y_range(1) & Y_vec <= y_range(2);

    x_fit = X_vec(in_region);
    y_fit = Y_vec(in_region);
    z_fit = Z_vec(in_region);

    % 删除 NaN 数据
    valid = ~isnan(z_fit);
    x_fit = x_fit(valid);
    y_fit = y_fit(valid);
    z_fit = z_fit(valid);

    if isempty(z_fit)
        fprintf('区域 %d 无有效数据，跳过。\n', k);
        continue;
    end

    % ====== 平面拟合：z = ax + by + c ======
    A = [x_fit, y_fit, ones(size(x_fit))];
    coeff = A \ z_fit;
    a = coeff(1); b = coeff(2); c = coeff(3);

    % ====== 坡度角计算 ======
    slope_angle = acos(1 / sqrt(a^2 + b^2 + 1));  % 弧度
    slope_angle_deg = rad2deg(slope_angle);

    % 输出信息
    fprintf('区域 %d: x ∈ [%.1f, %.1f], y ∈ [%.1f, %.1f]\n', ...
            k, x_range(1), x_range(2), y_range(1), y_range(2));
    fprintf('拟合平面：z = %.6f * x + %.6f * y + %.6f\n', a, b, c);
    fprintf('坡度角：%.5f°\n\n', slope_angle_deg);

    % ====== 可视化拟合平面 ======
    [xp, yp] = meshgrid(linspace(x_range(1), x_range(2), 20), ...
                        linspace(y_range(1), y_range(2), 20));
    zp = a * xp + b * yp + c;
    h = mesh(xp, yp, zp);
    set(h, 'FaceAlpha', 0.4);
    rectangle('Position', [x_range(1), y_range(1), ...
                           diff(x_range), diff(y_range)], ...
              'EdgeColor', 'r', 'LineWidth', 1);

    % ====== 测线布设 ======
    x_center = mean(x_range);              % 中心位置
    D = mean(z_fit) * 1852;                % 平均深度（米）
    D_nm = D / 1852;                       % 海里
    theta_deg = 120;
    target_overlap = 0.10;
    theta = deg2rad(theta_deg);
    phi1 = pi/2 - theta/2 - slope_angle;
    phi2 = pi/2 - theta/2 + slope_angle;

    % 测线间距 Δx 计算（起始间距）
    numerator_x1 = (D_nm * cos(slope_angle) * sin(theta / 2)) - 2 * sin(phi1);
    denominator_x1 = sin(phi1) - sin(slope_angle) * sin(theta / 2);
    delta_x = numerator_x1 / denominator_x1;

    % ====== 向右递推 ======
    x_right = x_center;
    x_right_list = x_right;

    for i = 1:100
        numerator_part1 = D_nm * (sin(phi1) + sin(phi2)) / (tan(slope_angle) * sin(phi2));
        coeff = sin(phi1) / (tan(slope_angle) * sin(theta / 2) * (1 - target_overlap));
        subtract = sin(phi1) / sin(phi2);
        numerator = numerator_part1 + x_right_list(end) * (coeff - subtract);
        denominator = 1 + coeff;
        x_next = numerator / denominator;

        if x_next > x_range(2)
            break;
        end
        x_right_list(end+1) = x_next;
    end

    % ====== 向左递推（对称） ======
    x_left_list = 2 * x_center - x_right_list(2:end);
    x_list = sort([x_left_list, x_right_list]);

    % ====== 输出测线位置 ======
    fprintf('区域 %d 共生成 %d 条测线（中心对称），单位：海里\n', k, length(x_list));
    disp(x_list);
end
