clear;
clc;
%% Grating Init
grt_l = 3;  % 相位突变点对数
grt_n = 2;  % 感兴趣的衍射级次, (2*grt_n + 1)对应目标等光强点数(分束比)
grt_spd = 1000;  % 采样点个数
grt_dc = 0.9;  % 预期总衍射效率

%% Optimizer Init
opt_tmax = 100;  % 初始温度
opt_tmin = 1.e-06;  % 终止温度
opt_decrate = 0.99;  % 下降速率
opt_its = 10 / grt_spd;  % 扰动强度
opt_alpha = 1.0;  % 误差函数的正则参数
opt_K = 1000;  % 迭代次数

%% Intermediate Variables Init
temp = opt_tmax;  % 当前迭代温度
a = sort(rand(1, 2*grt_l));  % 相位突变点坐标
t = disc_samp(a, grt_spd);  % 振幅透射率
p = diff_effi(t, grt_n);  % 各光强的衍射效率
eval = cost(p, opt_alpha, grt_dc);  % 代价
eval_all = zeros(1, 1000);  % 代价存储向量
count = 0;  % 计数器

%% Main
while(temp > opt_tmin)
    for i = 1:opt_K
        last_eval = eval;
        last_a = a;
        a = rand_dist(a, opt_its);
        t = disc_samp(a, grt_spd);
        p = diff_effi(t, grt_n);
        eval = cost(p, opt_alpha, grt_dc);
        if(eval < last_eval)
        elseif(exp(-(eval - last_eval) / temp) > 1 - 0.01 * rand())
        else
            a = last_a;
            eval = last_eval;
        end
    end
    temp = temp * opt_decrate;
    disp(eval);
    if(count < 999)
        count = count + 1;
        eval_all(count) = eval;
    end
end

%% Display
t = disc_samp(a, grt_spd);
x = linspace(0, 1, grt_spd);
figure;
plot(x, t);
title("振幅透射率系数");

figure;
plot(eval_all);
title("代价");

t_disp = disc_samp(a, 50);
T1 = repmat(t_disp, [50, 1]);
T = -T1 .* T1';
figure;
imshow(T, []);
title("图样");

T = repmat(T, [10, 10]);
T_f = fftshift(fft2(T));
T_fabs = abs(T_f) .^ 2;
T_fabs = T_fabs ./ (sum(T_fabs(:)));
figure;
imshow(T_fabs, []);
title("频谱");

%% For test

%% Discrete Sampling
% Input: 相位突变点坐标, (double) vector a; 采样点个数, (int) scalar spd
% Output: 密集采样后点列的振幅透射率系数, (double) vector t
function [t] = disc_samp(a, spd)
    t = zeros([1, spd]) - 1.;
    [~, col] = size(a);
    a = int32(a * spd);
    a(a <= 1) = 1;
    while(col)
        t(a(col-1):a(col)) = 1.;
        col = col - 2;
    end
end

%% Diffraction Efficiency Calculator
% Input: 密集采样后的振幅透射率系数, (double) vector t; 
%        感兴趣的衍射级次, (int) scalar n
% Output: 各级光强的衍射效率, (double) vector p, 对应于0,1,...,grt_n
function [p] = diff_effi(t, n) 
    t_f = fft(t);
    p_m = abs(t_f) .^ 2;
    p_m = p_m ./ sum(p_m);
    p = p_m(1:n+1);
end

%% Evaluation Function
% Input: 各级光强的衍射效率, (double) vector p; 正则参数, (double) scalar alpha;
%        预期总衍射效率, (double) scalar dc
% Output: 评价值, (double) scalar eval
function [eval] = cost(p, alpha, dc)
    [~, col] = size(p);
    pe = sum(p) + sum(p(2:end));
    p_hat = pe / (2 * col - 1);
    eval = (1 - alpha) * ((p(1) - p_hat) .^ 2 + 2 * sum((p(2:end) - p_hat) .^ 2)) + alpha * (1 - pe) .^ 2 ...
           + 10 * (sum(double(p > dc / (2 * col - 1))) + sum(double(p < 0.9 * dc / (2 * col - 1)))) * (p(1) - dc / (2 * col - 1)) .^ 2;
end

%% Random Disturbing
% Input: 相位突变点坐标, (double) vector a; 扰动强度, (double) scalar its
% Output: 扰动后的相位突变点坐标, (double) vector s
function [s] = rand_dist(a, its)
    [~, col] = size(a);
    a = (1 - 2 * rand(1,col)) .* its + a;
    count = 1;
    while(count < col)
        if(a(count + 1) < a(count))
            a(count + 1) = a(count);
        end
        count = count + 1;
    end
    a(a<0.) = 0.;
    a(a>1.) = 1.;
    s = sort(a);
end