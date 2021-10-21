close all; clear all;
oversample = 16;
Ne = 10;
figure
for i = 1:Ne
    S1 = decode(ones(1, 16));
    % tmp = decode(zeros(1, 18));
    % S2 = tmp(3:end);
    S2 = decode(zeros(1, 16));
    data = randi([0 1], [1, 32]);
    data(end-1:end) = [0, 0];
    DataPre = precode([S1, data, S2]);
    DataPre(47:48)
    CPM_BB = CPMmod(DataPre, oversample);
    plot(angle(CPM_BB))
    % plot(angle(CPM_BB))
    hold on;
    % plot(angle(pn_gen(S2)))
end

% data = randi([0 1], [1, 10]);
% % data = [0,1,1,0,1,0,0,1,1,1];
% data = [1,1,0,1,zeros(1,10)];
% CPM_BB = CPMmod(data, oversample);
% plot(angle(CPM_BB))