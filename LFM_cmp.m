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
CPM_BB = CPMmod(data, oversample);
GMSK_BB = GMSK_mode_new(data, oversample);

t = -0.5*Tp : 1/fs: 0.5*Tp - 1/fs;
t = t + 1/fs/2;
CPM_LFM = CPM_BB .* exp(1i*K*pi*t.^2);
GMSK_LFM = GMSK_BB .* exp(1i*K*pi*t.^2);

LFM = exp(1i*K*pi*t.^2);
LFM_AF = xcorr(LFM, LFM);
CPM_LFM_AF = xcorr(CPM_LFM, CPM_LFM);
GMSK_LFM_AF = xcorr(GMSK_LFM, GMSK_LFM);


figure
plot(20*log10(abs(LFM_AF)./max(abs(LFM_AF))));
hold on;
plot(20*log10(abs(CPM_LFM_AF)./max(abs(CPM_LFM_AF))));
hold on;
plot(20*log10(abs(GMSK_LFM_AF)./max(abs(GMSK_LFM_AF))));
ylim([-80 0]);
legend('LFM', 'CPM', 'GMSK')

Ne = 1000;

CPM_LFM_FFT = zeros(size(LFM));
GMSK_LFM_FFT = zeros(size(LFM));
for i = 1:Ne
    data = randi([0 1], [1 nn]);
    data(1:90) = 1;
    data(end-89:end) = 0;
    CPM_LFM = CPMmod(data, oversample) .* LFM;
    GMSK_LFM = GMSK_mode_new(data, oversample) .* LFM;
    CPM_LFM_FFT = CPM_LFM_FFT + abs(fftshift(fft(CPM_LFM)));
    GMSK_LFM_FFT = GMSK_LFM_FFT + abs(fftshift(fft(GMSK_LFM)));
end

CPM_LFM_FFT = CPM_LFM_FFT / Ne;
GMSK_LFM_FFT = GMSK_LFM_FFT / Ne;

% CPM_LFM_FFT = abs(fftshift(fft(CPM_LFM)));
% GMSK_LFM_FFT = abs(fftshift(fft(GMSK_LFM)));
LFM_FFT = abs(fftshift(fft(LFM)));

figure
plot(20*log10(CPM_LFM_FFT./max(CPM_LFM_FFT)))
hold on;
plot(20*log10(GMSK_LFM_FFT./max(GMSK_LFM_FFT)))
hold on;
plot(20*log10(LFM_FFT./max(LFM_FFT)))
ylim([-50 0]);
legend('CPM', 'GMSK', 'LFM')

