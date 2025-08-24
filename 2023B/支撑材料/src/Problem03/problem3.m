% 参数
D = 110;                    % 海水深度（米）
D_nm = D / 1852;
alpha_deg = 1.5;            % 坡度角（度）
theta_deg = 120;            % 开角（度）
target_overlap = 0.10;      % 目标重叠率（10%）
x_limit_nm = 2;             % 最大范围（海里）

% 弧度转换
alpha = deg2rad(alpha_deg);
theta = deg2rad(theta_deg);
phi1 = pi/2 - theta/2 - alpha;
phi2 = pi/2 - theta/2 + alpha;

% 初始化测线位置数组
x_list = [];

% 第一个测线位置（靠边界）
numerator_x1 = (D_nm * cos(alpha)* sin(theta / 2)) - 2 * sin(phi1);
denominator_x1 = sin(phi1) - sin(alpha) * sin(theta / 2);
x1 = numerator_x1 / denominator_x1;
x_list(1) = x1;

% 迭代计算后续测线位置

max_iter = 100;  % 防止死循环

i = 1;
while i <= max_iter
    % 当前测线位置
    xi = x_list(i);

    % 计算下一条测线位置
    % 计算中间项
    numerator_part1 = D_nm * (sin(phi1) + sin(phi2)) / (tan(alpha) * sin(phi2));

    coeff = sin(phi1) / (tan(alpha) * sin(theta / 2) * (1 - target_overlap));
    subtract = sin(phi1) / sin(phi2);

    numerator = numerator_part1 + xi * (coeff - subtract);
    denominator = 1 + coeff;

    % 递推下一条测线位置
    xi_next = numerator / denominator;


    % 判断是否超过最大范围
    if xi_next > x_limit_nm
        break;
    end

    % 存入结果
    x_list(i+1) = xi_next;
    i = i + 1;
end

% 输出
fprintf('共生成 %d 条测线，单位：海里\n', length(x_list));
disp('每条测线位置（单位：海里）:');
disp(x_list);

% ================== 绘制测线竖线图 ==================
figure;
hold on;

% 坐标范围
xlim([-2, 2]);
ylim([-1, 1]);

% 绘制每条测线（单位为海里，竖线）
for i = 1:length(x_list)
    x = x_list(i);  % 单位：海里
    line([x x], [-1, 1], 'Color', 'b', 'LineWidth', 1.5);  % 蓝色竖线
end

% 标注图形
xlabel('横向位置 x（海里）');
ylabel('纵向位置 y（海里）');
title('声呐测线布设示意图（竖线表示测线）');
grid on;
axis equal;

% 将测线位置加 2 海里偏移，并转换为米
x_offset_nm = 2;  % 海里偏移量
x_list_shifted_nm = x_list + x_offset_nm;     % 加偏移（单位：海里）
x_list_shifted_m = x_list_shifted_nm * 1852;  % 转换为米

% 输出结果
fprintf('\n每条测线位置（单位：米）:\n');
for i = 1:length(x_list_shifted_m)
    fprintf('第 %d 条测线位置：%.4f 米\n', i, x_list_shifted_m(i));
end

