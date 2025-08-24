[num,txt,raw] = xlsread('2021年国赛B题附件 (1).xlsx');

Y1 = num(:,2); % 第二列为乙醇转化率
Y2 = num(:,3); % 第三列为C4烯烃选择性
X4 = num(:,1); % 第一列为温度

% 初始化分组标记
group = ones(size(X4));
group_num = 1;
for i = 2:length(X4)
    if X4(i) < X4(i-1)
        group_num = group_num + 1;
    end
    group(i) = group_num;
end

num_groups = max(group);

for g = 1:num_groups
    if g <= 14
        group_label = ['A', num2str(g)];
    else
        group_label = ['B', num2str(g-14)];
    end
    
    idx = (group == g);
    x_group = X4(idx);
    y1_group = Y1(idx);
    y2_group = Y2(idx);

    % === 二次拟合 ===
    p1 = polyfit(x_group, y1_group, 2);
    y1_quad_fit = polyval(p1, x_group);
    R2_y1_quad = 1 - sum((y1_group - y1_quad_fit).^2)/sum((y1_group-mean(y1_group)).^2);

    p2 = polyfit(x_group, y2_group, 2);
    y2_quad_fit = polyval(p2, x_group);
    R2_y2_quad = 1 - sum((y2_group - y2_quad_fit).^2)/sum((y2_group-mean(y2_group)).^2);

    % === S 型拟合 ===
    logistic_fun = @(b,x) b(1) ./ (1 + exp(-b(2)*(x-b(3))));  % b = [L,k,x0]
    
    % 初始猜测参数
    beta0 = [max(y1_group), 0.1, mean(x_group)];
    b1 = nlinfit(x_group, y1_group, logistic_fun, beta0);
    y1_logi_fit = logistic_fun(b1, x_group);
    R2_y1_logi = 1 - sum((y1_group-y1_logi_fit).^2)/sum((y1_group-mean(y1_group)).^2);

    beta0 = [max(y2_group), 0.1, mean(x_group)];
    b2 = nlinfit(x_group, y2_group, logistic_fun, beta0);
    y2_logi_fit = logistic_fun(b2, x_group);
    R2_y2_logi = 1 - sum((y2_group-y2_logi_fit).^2)/sum((y2_group-mean(y2_group)).^2);

    % === 输出结果 ===
    fprintf('%s:\n', group_label);
    fprintf('  y1 二次 R^2 = %.4f, S型 R^2 = %.4f\n', R2_y1_quad, R2_y1_logi);
    fprintf('  y2 二次 R^2 = %.4f, S型 R^2 = %.4f\n', R2_y2_quad, R2_y2_logi);

    % === 画图 ===
    figure;
    hold on;
    x_fit = linspace(min(x_group), max(x_group), 100);

    % y1
    plot(x_group, y1_group, 'ro'); % 数据点
    plot(x_fit, polyval(p1, x_fit), 'r-', 'LineWidth',1.5); % 二次
    plot(x_fit, logistic_fun(b1,x_fit), 'r--', 'LineWidth',1.5); % Logistic

    % y2
    plot(x_group, y2_group, 'bo');
    plot(x_fit, polyval(p2, x_fit), 'b-', 'LineWidth',1.5);
    plot(x_fit, logistic_fun(b2,x_fit), 'b--', 'LineWidth',1.5);

    title([group_label,' 拟合对比']);
    xlabel('X4 (温度)');
    ylabel('Y');
    legend('y1 data','y1 二次','y1 S型','y2 data','y2 二次','y2 S型');
    grid on;
    hold off;
end
