function [h, w, t, t_str] = calculate_overlap_rate(D, y, alpha_deg, theta_deg)
    % 将输入的角度转换为弧度
    alpha = deg2rad(alpha_deg);
    theta = deg2rad(theta_deg);

    % 计算海水深度 h
    h = D - y * tan(alpha);

    % 计算覆盖宽度 w
    w = (D - y * tan(alpha)) * (1 / cos(alpha + theta / 2) + 1 / cos(alpha - theta / 2)) * sin(theta / 2);

    % 初始化重叠率数组
    t = NaN(size(w));

    % 计算重叠率 t
    for i = 2:length(w)
        % 计算相邻测线位置的差值
        delta_y = abs(y(i) - y(i-1));  % 使用绝对值计算相邻测线的距离

        term1 = h(i) * sin(theta / 2) / cos(alpha - theta / 2);
        term2 = (h(i) - delta_y * tan(alpha)) * sin(theta / 2) / cos(alpha + theta / 2);
        t(i) = (1 - delta_y / (term1 + term2)) * 100;
    end

    % 转换重叠率 NaN 为 '—'（用于文本显示）
    t_str = strings(size(t));
    for i = 1:length(t)
        if isnan(t(i))
            t_str(i) = '—';
        else
            t_str(i) = sprintf('%.4f', t(i));
        end
    end

    % 输出计算结果
    fprintf('测线距中心点处的距离/m:\n');
    disp(y);

    fprintf('海水深度/m:\n');
    disp(h);

    fprintf('覆盖宽度/m:\n');
    disp(w);

    fprintf('与前一条测线的重叠率/%%:\n');
    disp(t_str);
end
