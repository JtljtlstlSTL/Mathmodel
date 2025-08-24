clc;clear;
%  1. 读取数据 
data = xlsread('附件.xlsx');

% 提取坐标向量
x = data(1, 2:end);          % 横坐标向量
y = data(2:end, 1);          % 纵坐标向量
z_nm = data(2:end, 2:end);      % 深度矩阵
z = z_nm/1852;

% 构造原始网格
[X, Y] = meshgrid(x, y);

%  2. 插值生成更细网格 
[xi, yi] = meshgrid(linspace(min(x), max(x), 200), ...
                    linspace(min(y), max(y), 200));
zi = griddata(X, Y, z, xi, yi, 'cubic');

%  3. 绘制等深线图 
figure;
contourf(xi, yi, zi, 50);      % 填色等深线
colorbar;
hold on;
contour(xi, yi, zi, 50, 'k');  % 黑色轮廓线
xlabel('Longitude');
ylabel('Latitude');
title('等深线图及拟合平面叠加');
axis equal;

%  4. 指定拟合区域（手动指定） 
x_range = [0, 1.1];    % x 范围（可根据具体坐标调整）
y_range = [0, 2];    % y 范围

% 画出矩形边界
rectangle('Position', [x_range(1), y_range(1), ...
                       diff(x_range), diff(y_range)], ...
          'EdgeColor', 'r', 'LineWidth', 2);

%  5. 提取该区域内的插值点 
X_vec = xi(:);
Y_vec = yi(:);
Z_vec = zi(:);

in_region = X_vec >= x_range(1) & X_vec <= x_range(2) & ...
            Y_vec >= y_range(1) & Y_vec <= y_range(2);

x_fit = X_vec(in_region);
y_fit = Y_vec(in_region);
z_fit = Z_vec(in_region);

%  6. 最小二乘平面拟合 z = a*x + b*y + c 
A = [x_fit, y_fit, ones(size(x_fit))];
coeff = A \ z_fit;
a = coeff(1); b = coeff(2); c = coeff(3);

fprintf('拟合平面方程：z = %.4f * x + %.4f * y + %.4f\n', a, b, c);

%  7. 构造拟合平面网格 
[xp, yp] = meshgrid(linspace(x_range(1), x_range(2), 20), ...
                    linspace(y_range(1), y_range(2), 20));
zp = a * xp + b * yp + c;

%  8. 绘制拟合平面叠加 
mesh(xp, yp, zp);         % 绘制拟合平面
colormap winter
alpha(0.5)                % 设置透明度
