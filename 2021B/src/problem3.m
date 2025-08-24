% 定义目标函数（负 Y3，因为 fmincon 默认最小化）
objfun = @(X) -( (-80.526 -0.0337*X(1) +0.14124*X(2) -8.7659*X(3) +0.13327*X(4) +0.33321*X(5)) ...
                  * (-48.732 +0.0029892*X(1) +0.086447*X(2) +2.6734*X(3) -3.1639*X(4) +0.18122*X(5)) );


% 边界
lb = [33, 33, 0.5, 0.3, 250];
ub = [200, 200, 5, 2.1, 450];

% 初始猜测
x0 = mean([lb; ub]);

% 调用 fmincon
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x_opt, fval] = fmincon(objfun, x0, [], [], [], [], lb, ub, [], options);

% 输出最大 Y3
Y3_max = -fval;
Y3_pct = Y3_max / 100;
fprintf('最优自变量 X = [%f, %f, %f, %f, %f]\n', x_opt);
fprintf('最大 Y3 = %f\n', Y3_pct);