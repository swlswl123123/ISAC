close all; clear all;
open('CPM_LFM_40_sps.fig')
lh = findall(gca, 'type', 'line');
bit = get(lh, 'ydata');
x = get(lh, 'xdata');
close all;

ebn0 = x;
CPM_LFM = bit;
CPM_th = [0.0801170426,0.057440054,0.0383847914,0.0236406160,0.013104079,0.006416256,0.00271307059,0.000933942,0.0002597727,5.2137226e-5,7.1792964e-6];
MSK_th = [0.0786496035251426,0.0562819519765415,0.0375061283589260,0.0228784075610853,0.0125008180407376,0.00595386714777866,0.00238829078093281,0.000772674815378444,0.000190907774075993,3.36272284196176e-05,3.87210821552205e-06];
f = figure;
f.PaperUnits = 'centimeters';
f.PaperSize = [16, 12];
f.Units = 'centimeters';
f.Position = [0, 0, 16, 12];
semilogy(ebn0, CPM_LFM, '-^', 'LineWidth', 2);
hold on;
semilogy(ebn0, CPM_th, '-s', 'LineWidth', 2);
grid on;
% semilogy(ebn0, MSK_th, '-s', 'LineWidth', 2);
legend('k-LFM-CPM有频偏相偏','CPM理论性能')
set(get(gca, 'XLabel'), 'String', 'E_b/N_0');