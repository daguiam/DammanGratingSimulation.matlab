clear;
clc;
%% Grating Init
grt_l = 3;  % ��λͻ������
grt_n = 2;  % ����Ȥ�����伶��, (2*grt_n + 1)��ӦĿ��ȹ�ǿ����(������)
grt_spd = 1000;  % ���������
grt_dc = 0.9;  % Ԥ��������Ч��

%% Optimizer Init
opt_tmax = 100;  % ��ʼ�¶�
opt_tmin = 1.e-06;  % ��ֹ�¶�
opt_decrate = 0.99;  % �½�����
opt_its = 10 / grt_spd;  % �Ŷ�ǿ��
opt_alpha = 1.0;  % �������������
opt_K = 1000;  % ��������

%% Intermediate Variables Init
temp = opt_tmax;  % ��ǰ�����¶�
a = sort(rand(1, 2*grt_l));  % ��λͻ�������
t = disc_samp(a, grt_spd);  % ���͸����
p = diff_effi(t, grt_n);  % ����ǿ������Ч��
eval = cost(p, opt_alpha, grt_dc);  % ����
eval_all = zeros(1, 1000);  % ���۴洢����
count = 0;  % ������

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
title("���͸����ϵ��");

figure;
plot(eval_all);
title("����");

t_disp = disc_samp(a, 50);
T1 = repmat(t_disp, [50, 1]);
T = -T1 .* T1';
figure;
imshow(T, []);
title("ͼ��");

T = repmat(T, [10, 10]);
T_f = fftshift(fft2(T));
T_fabs = abs(T_f) .^ 2;
T_fabs = T_fabs ./ (sum(T_fabs(:)));
figure;
imshow(T_fabs, []);
title("Ƶ��");

%% For test

%% Discrete Sampling
% Input: ��λͻ�������, (double) vector a; ���������, (int) scalar spd
% Output: �ܼ���������е����͸����ϵ��, (double) vector t
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
% Input: �ܼ�����������͸����ϵ��, (double) vector t; 
%        ����Ȥ�����伶��, (int) scalar n
% Output: ������ǿ������Ч��, (double) vector p, ��Ӧ��0,1,...,grt_n
function [p] = diff_effi(t, n) 
    t_f = fft(t);
    p_m = abs(t_f) .^ 2;
    p_m = p_m ./ sum(p_m);
    p = p_m(1:n+1);
end

%% Evaluation Function
% Input: ������ǿ������Ч��, (double) vector p; �������, (double) scalar alpha;
%        Ԥ��������Ч��, (double) scalar dc
% Output: ����ֵ, (double) scalar eval
function [eval] = cost(p, alpha, dc)
    [~, col] = size(p);
    pe = sum(p) + sum(p(2:end));
    p_hat = pe / (2 * col - 1);
    eval = (1 - alpha) * ((p(1) - p_hat) .^ 2 + 2 * sum((p(2:end) - p_hat) .^ 2)) + alpha * (1 - pe) .^ 2 ...
           + 10 * (sum(double(p > dc / (2 * col - 1))) + sum(double(p < 0.9 * dc / (2 * col - 1)))) * (p(1) - dc / (2 * col - 1)) .^ 2;
end

%% Random Disturbing
% Input: ��λͻ�������, (double) vector a; �Ŷ�ǿ��, (double) scalar its
% Output: �Ŷ������λͻ�������, (double) vector s
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