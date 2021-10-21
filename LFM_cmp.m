clear all;close all;

% To Do: compare with MSK; calculate the rate of CPM(在频偏条件下完善解调); frequency reveal; CE-OFDM

% CPM parameter
h = 0.5; % modulation number
L = 2;   % related length
nn = 5000; % No. of modulate code
Tp = 50;  % LFM time width
Ts = Tp/nn;% symbol time width
sym_rate = 1/Ts;
% LFM parameter
B = 100; % Band width MHz
K = B/Tp;  % FM slope MHz/mus B/Tp
fc = 100; % carry frequence MHz 
oversample = 80;
fs = sym_rate * oversample; % sample frequence MHz 

data = randi([0 1], [1 nn]);
data(1:90) = 1;
data(end-89:end) = 0;

MSK_BB = MSKmod(data, oversample);
CPM_BB = CPMmod(data, oversample);
GMSK_BB = GMSK_mode_new(data, oversample);

t = -0.5*Tp : 1/fs: 0.5*Tp - 1/fs;
t = t + 1/fs/2;

LFM = exp(1i*K*pi*t.^2);
MSK_LFM = MSK_BB .* LFM;
CPM_LFM = CPM_BB .* LFM;
GMSK_LFM = GMSK_BB .* LFM;
LFM_AF = xcorr(LFM, LFM);

Ne = 1;

MSK_LFM_FFT = zeros(size(LFM));
CPM_LFM_FFT = zeros(size(LFM));
GMSK_LFM_FFT = zeros(size(LFM));

MSK_LFM_AF = zeros(size(LFM_AF));
CPM_LFM_AF = zeros(size(LFM_AF));
GMSK_LFM_AF = zeros(size(LFM_AF));

% signal = GMSK_LFM;
% signal_reverse = (fliplr(signal));
% % LFM_reverse_FFT = fft(LFM_reverse, 2*length(LFM)-1);
% signal_Ham = signal_reverse .* (hamming(length(signal)).');
% signal_AF = xcorr(signal, signal_Ham);

% figure
% plot(20*log10(abs(signal_AF)./max(abs(signal_AF))));
% hold on;

 
y=CPM_LFM;
yfft = fft(y, 2*Tp*fs-1);
h=zeros(1, Tp*fs);
for i=1:Tp*fs
    h(i)=conj(y(Tp*fs-i+1));
end
win = hamming(Tp*fs)';% Hamming 窗
hfft = fft(h .* win, 2*Tp*fs-1);
signal_AF = ifft(yfft .* hfft);
plot(20*log10(abs(signal_AF)./max(abs(signal_AF))));
hold on

% signal = LFM;
% signal_reverse = (fliplr(signal));
% % LFM_reverse_FFT = fft(LFM_reverse, 2*length(LFM)-1);
% signal_Ham = signal_reverse .* (hamming(length(signal)).');
% signal_AF = xcorr(signal, signal_Ham);
% plot(20*log10(abs(signal_AF)./max(abs(signal_AF))));



for i = 1:Ne
    data = randi([0 1], [1 nn]);
    data(1:90) = 1;
    data(end-89:end) = 0;
    MSK_LFM = MSKmod(data, oversample) .* LFM;
    CPM_LFM = CPMmod(data, oversample) .* LFM;
    GMSK_LFM = GMSK_mode_new(data, oversample) .* LFM;

    MSK_LFM_AF = MSK_LFM_AF + xcorr(MSK_LFM, MSK_LFM);
    CPM_LFM_AF = CPM_LFM_AF + xcorr(CPM_LFM, CPM_LFM);
    GMSK_LFM_AF = GMSK_LFM_AF + xcorr(GMSK_LFM, GMSK_LFM);

    MSK_LFM_FFT = MSK_LFM_FFT + abs(fftshift(fft(MSK_LFM)));
    CPM_LFM_FFT = CPM_LFM_FFT + abs(fftshift(fft(CPM_LFM)));
    GMSK_LFM_FFT = GMSK_LFM_FFT + abs(fftshift(fft(GMSK_LFM)));
end

MSK_LFM_AF = MSK_LFM_AF / Ne;
CPM_LFM_AF = CPM_LFM_AF / Ne;
GMSK_LFM_AF = GMSK_LFM_AF / Ne;

% figure
% plot(20*log10(abs(LFM_AF)./max(abs(LFM_AF))));
% hold on;
% plot(20*log10(abs(MSK_LFM_AF)./max(abs(MSK_LFM_AF))));
% hold on;
plot(20*log10(abs(CPM_LFM_AF)./max(abs(CPM_LFM_AF))));
% hold on;
% plot(20*log10(abs(GMSK_LFM_AF)./max(abs(GMSK_LFM_AF))));
% ylim([-80 0]);
% legend('LFM', 'MSK', 'CPM', 'GMSK')

% MSK_LFM_FFT = MSK_LFM_FFT / Ne;
% CPM_LFM_FFT = CPM_LFM_FFT / Ne;
% GMSK_LFM_FFT = GMSK_LFM_FFT / Ne;

% % CPM_LFM_FFT = abs(fftshift(fft(CPM_LFM)));
% % GMSK_LFM_FFT = abs(fftshift(fft(GMSK_LFM)));
% LFM_FFT = abs(fftshift(fft(LFM)));


% figure
% plot(20*log10(LFM_FFT./max(LFM_FFT)))
% hold on;
% plot(20*log10(MSK_LFM_FFT./max(MSK_LFM_FFT)))
% hold on;
% plot(20*log10(CPM_LFM_FFT./max(CPM_LFM_FFT)))
% hold on;
% plot(20*log10(GMSK_LFM_FFT./max(GMSK_LFM_FFT)))
% ylim([-50 0]);
% legend('LFM', 'MSK', 'CPM', 'GMSK')

