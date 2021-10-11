function [GMSK_BB] = GMSK_mode_new(data, oversamp)
%myFun - Description
%
% Syntax: [output] = myFun(input)
%
% Long description

bit_rate = 16e6;  % 符号速率
Tb = 1/bit_rate;  % 符号时间
BbTb = 0.3;       % 高斯滤波形状
Bb = BbTb / Tb;     
fs = bit_rate * oversamp; % 采样率
dt = 1/fs;        % 采样间隔

%% 生成高斯滤波函数gaussf与其积分表示
tt = -2.5*Tb : 1/fs: 2.5*Tb - 1/fs;
tt = tt + 1/fs/2;
gaussf = erfc(2 * pi * Bb * (tt - Tb / 2) / sqrt(log(2)) / sqrt(2)) / 2 - erfc(2 * pi * Bb * (tt + Tb / 2) / sqrt(log(2)) / sqrt(2)) / 2;
% plot(gaussf)
for i = 1:length(gaussf)

    if i == 1
        J_g(i) = gaussf(i) * dt;
    else
        J_g(i) = J_g(i - 1) + gaussf(i) * dt;
    end
end
J_g = J_g / 2 / Tb;

%% 调制GMSK
data = 2*data - 1;
bit_5 = [];
phi_all = [];
L = 0;
phi = zeros(1, oversamp);
for i = 1:length(data)
    if i == 1
        bit_5 = [0, 0, data(i:i+2)];
    elseif i == 2
        bit_5 = [0, data(i-1:i+2)];
    elseif i == length(data)-1
        bit_5 = [data(i-2:i+1), 0];
    elseif i == length(data)
        bit_5 = [data(i-2:i), 0, 0];
    else
        bit_5 = data(i-2:i+2);
    end
    if i >= 4
        L = L + data(i-3);
    end
    for j  = 1:oversamp
        phi(j) = bit_5(5)*J_g(j) + bit_5(4)*J_g(j + oversamp) + bit_5(3)*J_g(j + 2*oversamp) + bit_5(2)*J_g(j + 3*oversamp) + bit_5(1)*J_g(j + 4*oversamp);
    end
    phi_all = [phi_all, pi * phi + mod(L, 4) * pi / 2];
end
phi_all = mod(phi_all, 2*pi);
% t = 0:dt:length(data)*Tb-dt;
% plot(t, mod(phi_all, 2*pi)-2*pi)
GMSK_BB = exp(1i*phi_all);
end