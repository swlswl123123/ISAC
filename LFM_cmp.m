clear all;close all;

% To Do: compare with MSK; calculate the rate of CPM(在频偏条件下完善解调); frequency reveal; CE-OFDM

% CPM parameter
h = 0.5; % modulation number
L = 2;   % related length
nn = 500; % No. of modulate code
Tp = 50;  % LFM time width
Ts = Tp/nn;% symbol time width
sym_rate = 1/Ts;
% LFM parameter
B = 100; % Band width MHz
K = B/Tp;  % FM slope MHz/mus B/Tp
fc = 100; % carry frequence MHz 
oversample = 40;
fs = sym_rate * oversample; % sample frequence MHz 

data = randi([0 1], [1 nn]);
data(1:90) = 1;
data(end-89:end) = 0;

MSK_BB = MSKmod(data, oversample);
CPM_BB = CPMmod(data, oversample);
GMSK_BB = GMSK_mode_new(data, oversample);

t = -0.5*Tp : 1/fs: 0.5*Tp - 1/fs;
t = t + 1/fs/2;

MSK_LFM = MSK_BB .* exp(1i*K*pi*t.^2);
CPM_LFM = CPM_BB .* exp(1i*K*pi*t.^2);
GMSK_LFM = GMSK_BB .* exp(1i*K*pi*t.^2);

LFM = exp(1i*K*pi*t.^2);
LFM_AF = xcorr(LFM, LFM);
MSK_LFM_AF = xcorr(MSK_LFM, MSK_LFM);
CPM_LFM_AF = xcorr(CPM_LFM, CPM_LFM);
GMSK_LFM_AF = xcorr(GMSK_LFM, GMSK_LFM);

y=CPM_LFM;
yfft = fft(y, 2*Tp*fs-1);
% h=zeros(1, Tp*fs);
% for i=1:Tp*fs
%     h(i)=conj(y(Tp*fs-i+1));
% end
win = hamming(Tp*fs)';% Hamming 窗
winfft = fft(win, 2*Tp*fs-1);

CPM_LFM_AF = ifft(winfft.*yfft./yfft);

% hfft = fft(h.*win, 2*Tp*fs-1);     % 匹配滤波器的频域响应
% CPM_LFM_AF = abs(ifft(hfft.*yfft)); %脉冲压缩

% CPM_LFM_Reverse = fliplr(CPM_LFM);
% win = hamming(65);
% CPM_LFM_Reverse = conv(CPM_LFM_Reverse, win);
% CPM_LFM_Reverse = CPM_LFM_Reverse(32+1:end-32);
% CPM_LFM_AF = xcorr(CPM_LFM_Reverse, GMSK_LFM);



figure
% plot(20*log10(abs(LFM_AF)./max(abs(LFM_AF))));
% hold on;
% plot(20*log10(abs(MSK_LFM_AF)./max(abs(MSK_LFM_AF))));
% hold on;
plot(20*log10(abs(CPM_LFM_AF)./max(abs(CPM_LFM_AF))));
% hold on;
% plot(20*log10(abs(GMSK_LFM_AF)./max(abs(GMSK_LFM_AF))));
% ylim([-80 0]);
legend('LFM', 'MSK', 'CPM', 'GMSK')

Ne = 1;

MSK_LFM_FFT = zeros(size(LFM));
CPM_LFM_FFT = zeros(size(LFM));
GMSK_LFM_FFT = zeros(size(LFM));
for i = 1:Ne
    data = randi([0 1], [1 nn]);
    % data(1:90) = 1;
    % data(end-89:end) = 0;
    MSK_LFM = MSKmod(data, oversample) .* LFM;
    CPM_LFM = CPMmod(data, oversample) .* LFM;
    GMSK_LFM = GMSK_mode_new(data, oversample) .* LFM;
    MSK_LFM_FFT = MSK_LFM_FFT + abs(fftshift(fft(MSK_LFM)));
    CPM_LFM_FFT = CPM_LFM_FFT + abs(fftshift(fft(CPM_LFM)));
    GMSK_LFM_FFT = GMSK_LFM_FFT + abs(fftshift(fft(GMSK_LFM)));
end

MSK_LFM_FFT = MSK_LFM_FFT / Ne;
CPM_LFM_FFT = CPM_LFM_FFT / Ne;
GMSK_LFM_FFT = GMSK_LFM_FFT / Ne;

% CPM_LFM_FFT = abs(fftshift(fft(CPM_LFM)));
% GMSK_LFM_FFT = abs(fftshift(fft(GMSK_LFM)));
LFM_FFT = abs(fftshift(fft(LFM)));


figure
plot(20*log10(LFM_FFT./max(LFM_FFT)))
hold on;
plot(20*log10(MSK_LFM_FFT./max(MSK_LFM_FFT)))
hold on;
plot(20*log10(CPM_LFM_FFT./max(CPM_LFM_FFT)))
hold on;
plot(20*log10(GMSK_LFM_FFT./max(GMSK_LFM_FFT)))
ylim([-50 0]);
legend('LFM', 'MSK', 'CPM', 'GMSK')

