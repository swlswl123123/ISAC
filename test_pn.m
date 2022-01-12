close all; clear all;
oversample = 128;
Ne = 100;
figure
x = 0:1/128:500 - 1/128;
for i = 1:Ne

    S1 = ones(1, 16);
    S1(2:2:end) = 0;
    % tmp = decode(zeros(1, 18));
    % S2 = tmp(3:end);
    S2 = S1;

    % S1 = decode(ones(1, 16));
    % S2 = decode(zeros(1, 16));


    data = randi([0 1], [1, 130]);
    data(end-1:end) = [1, 0];
    DataPre = precode([S1, data, S2]);
%     DataPre(47:48)
    CPM_BB = CPMmod(DataPre, oversample);
    LFM = exp(1i*2*pi*x.^2);
    plot(angle(CPM_BB))
    hold on;
    % figure
    % subplot(2,1,1)
    % plot(x(1:32*oversample), real(CPM_BB(1:32*oversample).*LFM(1:32*oversample)))
    % xlim([1 32])
    % set(get(gca, 'Title'), 'String', 'I路');
    % subplot(2,1,2)
    % plot(x(1:32*oversample), imag(CPM_BB(1:32*oversample).*LFM(1:32*oversample)))
    % xlim([1 32])
    % set(get(gca, 'Title'), 'String', 'Q路');
    % plot(angle(pn_gen(S2)))
end

% data = randi([0 1], [1, 10]);
% % data = [0,1,1,0,1,0,0,1,1,1];
% data = [1,1,0,1,zeros(1,10)];
% CPM_BB = CPMmod(data, oversample);
% plot(angle(CPM_BB))
% xlabel('符号个数')
% ylabel('相位')