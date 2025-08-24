clc; clear;


data = xlsread('附件.xlsx');

% 提取坐标向量
x = data(1, 2:end);              % 横坐标（单位：海里）
y = data(2:end, 1);              % 纵坐标（单位：海里）
z_nm = data(2:end, 2:end);       % 深度（单位：米）
z = z_nm / 1852;                 % 转换为海里

% 构造原始网格
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

    % 平面拟合：z = ax + by + c
    A = [x_fit, y_fit, ones(size(x_fit))];
    coeff = A \ z_fit;  % 求解最小二乘法的系数
    a = coeff(1); b = coeff(2); c = coeff(3);

    % 预测值（拟合面上的 z）
    z_pred = A * coeff;  % 根据拟合系数得到的拟合平面值

    % 残差：实际值与拟合值之间的差距
    residual = z_fit - z_pred;

    % 计算均方误差（MSE）
    mse = mean(residual.^2);  % 平均平方误差
    rmse = sqrt(mse);         % 均方根误差

    % 计算决定系数 R²
    z_mean = mean(z_fit);  % 实际值的平均值
    ss_total = sum((z_fit - z_mean).^2);  % 总平方和
    ss_res = sum(residual.^2);  % 残差平方和
    R2 = 1 - ss_res / ss_total;  % 决定系数 R²

    % 计算坡度角（cos α = 1 / sqrt(a^2 + b^2 + 1)）
    slope_angle = acos(1 / sqrt(a^2 + b^2 + 1));  % 坡度角（弧度）
    slope_angle_deg = rad2deg(slope_angle);       % 转为度数

    % 输出拟合信息
    fprintf('区域 %d: x ∈ [%.1f, %.1f], y ∈ [%.1f, %.1f]\n', ...
            k, x_range(1), x_range(2), y_range(1), y_range(2));
    fprintf('拟合平面：z = %.6f * x + %.6f * y + %.6f\n', a, b, c);
    fprintf('拟合精度：RMSE = %.5f，R² = %.5f\n', rmse, R2);
    fprintf('坡度角：%.5f°\n\n', slope_angle_deg);

    % 可视化拟合平面
    [xp, yp] = meshgrid(linspace(x_range(1), x_range(2), 20), ...
                        linspace(y_range(1), y_range(2), 20));
    zp = a * xp + b * yp + c;
    
    h = mesh(xp, yp, zp);         % 绘制平面网格
    set(h, 'FaceAlpha', 0.4);     % 设置透明度
    rectangle('Position', [x_range(1), y_range(1), ...
                           diff(x_range), diff(y_range)], ...
              'EdgeColor', 'r', 'LineWidth', 1);
end
