% 读取数据
data = xlsread('附件.xlsx');

% 提取坐标向量
x = data(1, 2:end);          % 第一行，第2列到最后列，是x坐标（横向）
y = data(2:end, 1);          % 第二行到最后行，第1列，是y坐标（纵向）

% 提取深度矩阵 z
z = data(2:end, 2:end);      % 去掉第一行和第一列的纯z值矩阵

% 构造原始网格
[X, Y] = meshgrid(x, y);

% 目标插值网格（更细致）
[xi, yi] = meshgrid(linspace(min(x), max(x), 200), linspace(min(y), max(y), 200));

% 插值（注意用 meshgrid 格式）
zi = griddata(X, Y, z, xi, yi, 'cubic');

% 绘图
figure;
contourf(xi, yi, zi, 50);      % 填色等高线
colorbar;
hold on;
contour(xi, yi, zi, 50, 'k');  % 黑色轮廓线
xlabel('Longitude');
ylabel('Latitude');
title('等深线图（Bathymetric Contour Map）');
axis equal;
