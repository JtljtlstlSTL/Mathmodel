[num,txt,raw] = xlsread('2021年国赛B题附件 (1).xlsx');

Y1 = num(:,2); % 第二列为乙醇转化率
Y2 = num(:,3); % 第四列为C4烯烃选择性

X4 = num(:,1); % 第一列为温度

% 边界修正
Y1(Y1 <= 0) = 0.0001;
Y1(Y1 >= 100) = 99.9999;
Y2(Y2 <= 0) = 0.0001;
Y2(Y2 >= 100) = 99.9999;

% 计算新变量
y1 = log(Y1 ./ (100 - Y1));
y2 = log(Y2 ./ (100 - Y2));
y3 = log(sqrt(Y1 .* Y2) ./ (100 - sqrt(Y1 .* Y2)));

result = [X4, y1, y2];

% disp(result);


% delete('problem1_logit.xlsx');  % 删除旧文件
% xlswrite('problem1_logit.xlsx', result, 'Sheet1','A1');






% 初始化分组标记
group = ones(size(X4)); % 第一行属于第一组
group_num = 1;

for i = 2:length(X4)
    if X4(i) < X4(i-1)  % 温度下降，分组号加1
        group_num = group_num + 1;
    end
    group(i) = group_num;
end

% 合并结果
result_with_group = [result, group];

% disp(result_with_group);





num_groups = max(group);
coeff_all = zeros(num_groups, 2);

for g = 1:num_groups

    % 生成组标签
    if g <= 14
        group_label = ['A', num2str(g)];
    else
        group_label = ['B', num2str(g-14)];
    end

    idx = (group == g); % 当前组索引
    x_group = X4(idx);
    
    %---------------- y1 拟合 ----------------
    y_group1 = y1(idx);
    p1 = polyfit(x_group, y_group1, 2);
    coeff(g).y1 = p1;
    y1_fit_group = polyval(p1, x_group);

    % R²计算
    SS_res1 = sum((y_group1 - y1_fit_group).^2);
    SS_tot1 = sum((y_group1 - mean(y_group1)).^2);
    R2_y1 = 1 - SS_res1/SS_tot1;

    %---------------- y2 拟合 ----------------
    y_group2 = y2(idx);
    p2 = polyfit(x_group, y_group2, 2);
    coeff(g).y2 = p2;
    y2_fit_group = polyval(p2, x_group);

    % R²计算
    SS_res2 = sum((y_group2 - y2_fit_group).^2);
    SS_tot2 = sum((y_group2 - mean(y_group2)).^2);
    R2_y2 = 1 - SS_res2/SS_tot2;

    %---------------- 绘图 ----------------
    figure;
    hold on;
    x_fit = linspace(min(x_group), max(x_group), 50);  % 拟合曲线点

    % y1 曲线
    y1_fit_plot = polyval(p1, x_fit);
    scatter(x_group, y_group1, 36, 'r', 'filled');   % 原始数据
    plot(x_fit, y1_fit_plot, 'r-', 'LineWidth', 1.5);    % 拟合曲线

    % y2 曲线
    y2_fit_plot = polyval(p2, x_fit);
    scatter(x_group, y_group2, 36, 'b', 'filled');   % 原始数据
    plot(x_fit, y2_fit_plot, 'b-', 'LineWidth', 1.5);    % 拟合曲线

    title([group_label, ' 二次拟合 y1 & y2']);
    xlabel('X4 (温度)');
    ylabel('y');
    legend('y1 data','y1 fit','y2 data','y2 fit');
    grid on;
    hold off;

    % Pearson 相关性
    [R1,P1] = corrcoef(x_group, y_group1);
    [R2,P2] = corrcoef(x_group, y_group2);

    % 输出结果
    fprintf('%s:\n', group_label);
    fprintf('  y1 vs X4: r = %.4f, p = %.4f, R² = %.4f\n', R1(1,2), P1(1,2), R2_y1);
    fprintf('  y2 vs X4: r = %.4f, p = %.4f, R² = %.4f\n', R2(1,2), P2(1,2), R2_y2);

    %---------------- 绘图（反变换回原始 Y 空间） ----------------
    figure;
    hold on;
    x_fit = linspace(min(x_group), max(x_group), 100);

    % y1 曲线（反logit）
    y1_fit_plot = polyval(p1, x_fit);
    Y1_fit_plot = 100 ./ (1 + exp(-y1_fit_plot));

    scatter(x_group, Y1(idx), 36, 'r', 'filled'); % 原始 Y1 数据
    plot(x_fit, Y1_fit_plot, 'r-', 'LineWidth', 1.5); % 拟合曲线

    % y2 曲线（反logit）
    y2_fit_plot = polyval(p2, x_fit);
    Y2_fit_plot = 100 ./ (1 + exp(-y2_fit_plot));

    scatter(x_group, Y2(idx), 36, 'b', 'filled'); % 原始 Y2 数据
    plot(x_fit, Y2_fit_plot, 'b-', 'LineWidth', 1.5); % 拟合曲线

    title([group_label, ' Logistic反变换拟合 Y1 & Y2']);
    xlabel('X4 (温度)');
    ylabel('Y (%)');
    legend('Y1 data','Y1 fit','Y2 data','Y2 fit');
    grid on;
    hold off;

end

for g = 1:num_groups

     % 生成组标签
    if g <= 14
        group_label = ['A', num2str(g)];
    else
        group_label = ['B', num2str(g-14)];
    end
    disp([group_label, ' y1系数: ', num2str(coeff(g).y1)]);
    disp([group_label, ' y2系数: ', num2str(coeff(g).y2)]);
end



