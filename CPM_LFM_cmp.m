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
oversample = 8;
fs = sym_rate * oversample; % sample frequence MHz 

t = -0.5*Tp : 1/fs: 0.5*Tp - 1/fs;
t = t + 1/fs/2;

LFM = exp(1i*K*pi*t.^2);

Ne = 1;

CPM_LFM_Ham_AF = zeros(1, 2*Tp*fs-1);
CPM_LFM_AF = zeros(1, 2*Tp*fs-1);
CPM_LFM_FFT = zeros(1, Tp*fs);
CPM_LFM_origin_FFT = zeros(1, Tp*fs);

for i = 1:Ne
    data = randi([0 1], [1 nn]);
    data_origin = data;
    data(1:180) = 1;
    data(end-179:end) = 0;
    CPM_LFM = CPMmod(data, oversample) .* LFM;
    CPM_LFM_origin = CPMmod(data_origin, oversample) .* LFM;

    CPM_LFM_AF = CPM_LFM_AF + xcorr(CPM_LFM, CPM_LFM);
    y=CPM_LFM;
    yfft = fft(y, 2*Tp*fs-1);
    h=zeros(1, Tp*fs);
    for i=1:Tp*fs
        h(i)=conj(y(Tp*fs-i+1));
    end
    win = hamming(Tp*fs)';% Hamming 窗
    hfft = fft(h .* win, 2*Tp*fs-1);
    CPM_LFM_Ham_AF = CPM_LFM_Ham_AF + ifft(yfft .* hfft);;

    CPM_LFM_FFT = CPM_LFM_FFT + abs(fftshift(fft(CPM_LFM)));
    CPM_LFM_origin_FFT = CPM_LFM_origin_FFT + abs(fftshift(fft(CPM_LFM_origin)));
end

CPM_LFM_AF = CPM_LFM_AF / Ne;
CPM_LFM_Ham_AF = CPM_LFM_Ham_AF / Ne;

figure
plot(20*log10(abs(CPM_LFM_AF)./max(abs(CPM_LFM_AF))));
hold on;
plot(20*log10(abs(CPM_LFM_Ham_AF)./max(abs(CPM_LFM_Ham_AF))));

figure
plot(20*log10(abs(CPM_LFM_FFT)./max(abs(CPM_LFM_FFT))));
hold on;
plot(20*log10(abs(CPM_LFM_origin_FFT)./max(abs(CPM_LFM_origin_FFT))));