clear all; close all;
fs = 320;
B = 40;
T = 20;
k = B / T;
t = -0.5*T:1/fs:0.5*T-1/fs;

TotalNum = 800;
HeadLen = 128;
Ts = T / TotalNum; % 0.025
oversample = fs*Ts; % 8
lfm_phase = pi*k*t.^2;

Data = randi([0 1], [1 TotalNum]);
S1 = ones(1,HeadLen);
S1(2:2:end) = 0;
Data(1:HeadLen) = S1;
Data(end-HeadLen+1:end) = S1;

CPM_BB = CPMmod(Data, oversample);
CPM_LFM = CPM_BB .* exp(1i*lfm_phase);

wavHead = CPM_LFM(1:HeadLen*oversample);
signal_recv = [zeros(1, 1500), CPM_LFM, zeros(1, 1500)];

cor = zeros(1,3000);
for i = 1:length(signal_recv)-HeadLen*oversample
    cor(i) = abs(signal_recv(i:i+HeadLen*oversample-1) * wavHead');
end

plot(cor)

% cpm_phase = angle(CPMmod(Data, oversample));
% lfm_cpm_phase = lfm_phase + cpm_phase;
% I = cos(lfm_cpm_phase);
% Q = sin(lfm_cpm_phase);
