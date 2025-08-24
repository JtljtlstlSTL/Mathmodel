[num,txt,raw] = xlsread('2021年国赛B题附件 (2).xlsx');

Y = num(:,3:end);        % 选择性数据 (10行 × 6列)

% 转置：6个选择性 → 样本，10个时间点 → 特征
Y_t = Y';   % (6行 × 10列)

% 标准化（按每种选择性）
Y_norm = normalize(Y_t, 2);

% K-means 聚类（分2类）
k = 2;
[idx,C] = kmeans(Y_norm, k, 'Replicates',20);

% 获取选择性名称
sel_names = txt(3,3:end)';  % 假设第一行是表头，列3开始是选择性名称

% 输出结果
disp('各选择性所属的类别：');
for i = 1:size(Y_t,1)
    fprintf('%s -> 类别 %d\n', sel_names{i}, idx(i)); 
end


% 定义实际时间点
time_points = [20 70 110 163 197 240 273];

% 可视化：聚类结果热力图（类别1）
figure;
imagesc(Y_norm(idx==1,:));
colorbar;
xticks(1:length(time_points));
xticklabels(time_points);      % 时间点
yticks(1:sum(idx==1));
yticklabels(sel_names(idx==1));     % 选择性名字
xlabel('时间点 (min)');
ylabel('选择性');
title('类别1 (归一化后的选择性随时间变化)');

% 类别2
figure;
imagesc(Y_norm(idx==2,:));
colorbar;
xticks(1:length(time_points));
xticklabels(time_points);      % 时间点
yticks(1:sum(idx==2));
yticklabels(sel_names(idx==2));     % 选择性名字
xlabel('时间点 (min)');
ylabel('选择性');
title('类别2 (归一化后的选择性随时间变化)');
