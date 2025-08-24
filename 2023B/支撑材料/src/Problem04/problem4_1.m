% 读取数据
data = xlsread('附件.xlsx');

% 提取坐标向量
x = data(1, 2:end);
y = data(2:end, 1);
z = data(2:end, 2:end);

% 构造原始网格
[X, Y] = meshgrid(x, y);

% 插值网格
[xi, yi] = meshgrid(linspace(min(x), max(x), 200), ...
                    linspace(min(y), max(y), 200));
zi = griddata(X, Y, z, xi, yi, 'cubic');

% 绘图：只画等深线，不填色
figure;
contour(xi, yi, zi, 50, 'k');  % 黑色等深线
xlabel('Longitude');
ylabel('Latitude');
title('等深线图');
axis equal;
grid on;
