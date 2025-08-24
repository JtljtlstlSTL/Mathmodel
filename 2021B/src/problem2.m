% X1为Co/SiO2 质量，X2为HAP 质量，X3为乙醇浓度，X4为Co 负载量，X5为温度

%% 读取数据
[num, txt, raw] = xlsread('2021年国赛B题附件 (1).xlsx');

X5 = num(:,1);        % 温度
Y1 = num(:,2);        % 乙醇转化率
Y2 = num(:,4);        % C4烯烃选择性
catalyst_raw = raw(:,2);  % 催化剂组合在第2列

%% Step 1: 按温度递减分组
group = ones(size(X5));
group_raw_row = 2;
group_num = 1;
for i = 2:length(X5)
    if X5(i) < X5(i-1)
        group_num = group_num + 1;
        group_raw_row = [group_raw_row, i+1];
    end
    group(i) = group_num;
end
numGroups = group_num;



%% Step 2: 提取催化剂组合，并解析为 X1-X4
X1 = zeros(size(X5));  % Co/SiO2 质量 mg
X2 = zeros(size(X5));  % HAP 质量 mg
X3 = zeros(size(X5));  % 乙醇流量 ml/min
X4 = zeros(size(X5));  % Co 负载量 wt%

for i = 1:length(X4)
    j = group_raw_row(group(i));
    text = catalyst_raw{j};  % 例如 '200mg 1wt%Co/SiO2-200mg HAP-乙醇浓度1.68ml/min'
    
    % 提取 Co/SiO2 质量
    tok = regexp(text, '(\d+\.?\d*)mg \d+\.?\d*wt%Co/SiO2', 'tokens');
    if ~isempty(tok), X1(i) = str2double(tok{1}{1}); end
    
    % 提取 HAP 质量
    tok = regexp(text, '(\d+\.?\d*)mg HAP', 'tokens');
    if ~isempty(tok)
        X2(i) = str2double(tok{1}{1});
    else
        X2(i) = 0;  % 无HAP
    end
    
    % 提取乙醇流量
    tok = regexp(text, '乙醇浓度(\d+\.?\d*)ml/min', 'tokens');
    if ~isempty(tok), X3(i) = str2double(tok{1}{1}); end
    
    % 提取 Co 负载量
    tok = regexp(text, '(\d+\.?\d*)wt%Co/SiO2', 'tokens');
    if ~isempty(tok), X4(i) = str2double(tok{1}{1}); end
end

%% Step 3: 构建设计矩阵
X = [X1, X2, X3, X4, X5];

%% Step 4: 多元回归
mdl1 = fitlm(X, Y1, 'linear');
mdl2 = fitlm(X, Y2, 'linear');

%% Step 5: 可视化
figure;
subplot(1,2,1); plot(mdl1); title('乙醇转化率 多元回归');
subplot(1,2,2); plot(mdl2); title('C4烯烃选择性 多元回归');

%% Step 6: 导出回归系数表格
tbl1 = mdl1.Coefficients;  % 包含 Estimate, SE, tStat, pValue
tbl2 = mdl2.Coefficients;

writetable(tbl1, '回归结果_乙醇转化率.xlsx');
writetable(tbl2, '回归结果_C4烯烃选择性.xlsx');

%% Step 7: 打印回归系数分析表
disp('===== 乙醇转化率 回归系数分析 =====');
disp(tbl1);

disp('===== C4烯烃选择性 回归系数分析 =====');
disp(tbl2);

%% Step 8: 输出模型指标
fprintf('乙醇转化率: R²=%.4f, 调整R²=%.4f, F=%.4f, p=%.4f\n', ...
    mdl1.Rsquared.Ordinary, mdl1.Rsquared.Adjusted, ...
    mdl1.ModelFitVsNullModel.Fstat, mdl1.ModelFitVsNullModel.Pvalue);

fprintf('C4烯烃选择性: R²=%.4f, 调整R²=%.4f, F=%.4f, p=%.4f\n', ...
    mdl2.Rsquared.Ordinary, mdl2.Rsquared.Adjusted, ...
    mdl2.ModelFitVsNullModel.Fstat, mdl2.ModelFitVsNullModel.Pvalue);