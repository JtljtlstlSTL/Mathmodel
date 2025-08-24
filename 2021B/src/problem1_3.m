[num,txt,raw] = xlsread('2021年国赛B题附件 (2).xlsx');

Y = num(:,2:end);        % 选择性数据 (10行 × 6列)
sel_names = [{txt{2,2}}; txt(3,3:end)'];  
time_points = [20 70 110 163 197 240 273];  % 实际时间点

figure; hold on;
colors = lines(size(Y,2));

for i = 1:size(Y,2)
    y_orig = Y(:,i);  % 原始百分比选择性
    x = time_points';
    
    % Logit 变换
    y_logit = log(y_orig ./ (100 - y_orig));
    
    % 二次拟合
    p = polyfit(x, y_logit, 2);
    
    % 拟合值（logit）
    y_logit_fit = polyval(p, x);
    
    % 反变换回百分比
    y_fit = 100 ./ (1 + exp(-y_logit_fit));
    
    % === 计算 R² ===
    SS_res = sum((y_orig - y_fit).^2);
    SS_tot = sum((y_orig - mean(y_orig)).^2);
    R2 = 1 - SS_res/SS_tot;
    
    % 输出拟合方程 + R²
    fprintf('%s 拟合方程 (Logit二次)：\n', sel_names{i});
    fprintf('logit(Y) = %.8f * x^2 + %.8f * x + %.8f\n', p(1), p(2), p(3));
    fprintf('反变换百分比形式：Y = 100 / (1 + exp(-(%.8f * x^2 + %.8f * x + %.8f)))\n', p(1), p(2), p(3));
    fprintf('R² = %.4f\n\n', R2);
    
    % 拟合点（细分时间）
    x_fine = linspace(min(x), max(x), 200);
    y_fit_logit = polyval(p, x_fine);
    y_fit_curve = 100 ./ (1 + exp(-y_fit_logit));
    
    % 绘制原始折线
    plot(x, y_orig, 'o-', 'Color', colors(i,:), 'DisplayName', [sel_names{i} ' 原始']);
    
    % 绘制拟合曲线
    plot(x_fine, y_fit_curve, '--', 'Color', colors(i,:), 'DisplayName', [sel_names{i} ' 拟合']);
end

xlabel('时间 (min)');
ylabel('选择性 (%)');
title('各选择性原始数据折线 + Logit 二次拟合曲线');
legend('Location','best');
grid on;
hold off;
