function out = CPMdemodD2(CPM_recv)
%CPMdemod - Description
%
% Syntax: out = CPMdemod(CPM_recv)
%
% Long description
nn = 50;
path_rec = cell(8, 1);
weight = zeros(8, nn);

down_sample = 4;

phi_000_000 = exp(-1i * pi) * ones(1, down_sample);
phi_000_001 = exp(1i * (((1:down_sample)-1) * (pi/2/down_sample) - pi));
phi_001_010 = exp(1i * (((1:down_sample)-1) * (pi/2/down_sample) - pi/2));
phi_001_011 = exp(1i * (((1:down_sample)-1) * (pi/down_sample) - pi/2));
phi_010_100 = exp(-1i * ((1:down_sample)-1) * (pi/2/down_sample));
phi_010_101 = ones(1, down_sample);
phi_011_110 = exp(1i * pi/2) * ones(1, down_sample);
phi_011_111 = exp(1i* (((1:down_sample)-1) * (pi/2/down_sample) + pi/2));
phi_100_000 = exp(1i* (-((1:down_sample)-1) * (pi/2/down_sample) - pi/2));
phi_100_001 = exp(-1i * pi/2) * ones(1, down_sample);
phi_101_010 = ones(1, down_sample);
phi_101_011 = exp(1i * ((1:down_sample)-1) * (pi/2/down_sample));
phi_110_100 = exp(1i * (-((1:down_sample)-1) * (pi/down_sample) + pi/2));
phi_110_101 = exp(1i * (-((1:down_sample)-1) * (pi/2/down_sample) + pi/2));
phi_111_110 = exp(1i * (-((1:down_sample)-1) * (pi/2/down_sample) + pi));
phi_111_111 = exp(1i * pi) * ones(1, down_sample);

for i = 3:nn
    CPM_dif = CPM_recv((i-1)*down_sample+1:i*down_sample) .* conj(CPM_recv((i-3)*down_sample+1:(i-2)*down_sample));
    if i == 3
        phi_000_s = exp(1i * (-((1:down_sample)-1) * (pi/4/down_sample) - 3*pi/4));
        phi_001_s = exp(1i * (((1:down_sample)-1) * (pi/4/down_sample) - 3*pi/4));
        phi_010_s = exp(1i * (((1:down_sample)-1) * (pi/4/down_sample) - pi/4));
        phi_011_s = exp(1i * (((1:down_sample)-1) * (3*pi/4/down_sample) - pi/4));
        phi_100_s = exp(1i * (-((1:down_sample)-1) * (3*pi/4/down_sample) + pi/4));
        phi_101_s = exp(1i * (-((1:down_sample)-1) * (pi/4/down_sample) + pi/4));
        phi_110_s = exp(1i * (-((1:down_sample)-1) * (pi/4/down_sample) + 3*pi/4));
        phi_111_s = exp(1i * (((1:down_sample)-1) * (pi/4/down_sample) + 3*pi/4));
        weight(1, i) = real(sum(CPM_dif.*conj(phi_000_s)));
        weight(2, i) = real(sum(CPM_dif.*conj(phi_001_s)));
        weight(3, i) = real(sum(CPM_dif.*conj(phi_010_s)));
        weight(4, i) = real(sum(CPM_dif.*conj(phi_011_s)));
        weight(5, i) = real(sum(CPM_dif.*conj(phi_100_s)));
        weight(6, i) = real(sum(CPM_dif.*conj(phi_101_s)));
        weight(7, i) = real(sum(CPM_dif.*conj(phi_110_s)));
        weight(8, i) = real(sum(CPM_dif.*conj(phi_111_s)));
        path_rec{1} = [0, 0, 0];
        path_rec{2} = [0, 0, 1];
        path_rec{3} = [0, 1, 0];
        path_rec{4} = [0, 1, 1];
        path_rec{5} = [1, 0, 0];
        path_rec{6} = [1, 0, 1];
        path_rec{7} = [1, 1, 0];
        path_rec{8} = [1, 1, 1];
    else
        rec1 = [];
        rec2 = [];
        rec3 = [];
        rec4 = [];
        rec5 = [];
        rec6 = [];
        rec7 = [];
        rec8 = [];

        w_000_000 = weight(1, i-1) + real(sum(CPM_dif.*conj(phi_000_000)));
        w_000_001 = weight(1, i-1) + real(sum(CPM_dif.*conj(phi_000_001)));
        w_001_010 = weight(2, i-1) + real(sum(CPM_dif.*conj(phi_001_010)));
        w_001_011 = weight(2, i-1) + real(sum(CPM_dif.*conj(phi_001_011)));
        w_010_100 = weight(3, i-1) + real(sum(CPM_dif.*conj(phi_010_100)));
        w_010_101 = weight(3, i-1) + real(sum(CPM_dif.*conj(phi_010_101)));
        w_011_110 = weight(4, i-1) + real(sum(CPM_dif.*conj(phi_011_110)));
        w_011_111 = weight(4, i-1) + real(sum(CPM_dif.*conj(phi_011_111)));

        w_100_000 = weight(5, i-1) + real(sum(CPM_dif.*conj(phi_100_000)));
        w_100_001 = weight(5, i-1) + real(sum(CPM_dif.*conj(phi_100_001)));
        w_101_010 = weight(6, i-1) + real(sum(CPM_dif.*conj(phi_101_010)));
        w_101_011 = weight(6, i-1) + real(sum(CPM_dif.*conj(phi_101_011)));
        w_110_100 = weight(7, i-1) + real(sum(CPM_dif.*conj(phi_110_100)));
        w_110_101 = weight(7, i-1) + real(sum(CPM_dif.*conj(phi_110_101)));
        w_111_110 = weight(8, i-1) + real(sum(CPM_dif.*conj(phi_111_110)));
        w_111_111 = weight(8, i-1) + real(sum(CPM_dif.*conj(phi_111_111)));

        if w_000_000 > w_100_000
            weight(1, i) = w_000_000;
            rec1 = [path_rec{1}, 0];
        else
            weight(1, i) = w_100_000;
            rec1 = [path_rec{5}, 0];
        end

        if w_000_001 > w_100_001
            weight(2, i) = w_000_001;
            rec2 = [path_rec{1}, 1];
        else
            weight(2, i) = w_100_000;
            rec2 = [path_rec{5}, 1];
        end

        if w_001_010 > w_101_010
            weight(3, i) = w_001_010;
            rec3 = [path_rec{2}, 0];
        else
            weight(3, i) = w_101_010;
            rec3 = [path_rec{6}, 0];
        end

        if w_001_011 > w_101_011
            weight(4, i) = w_001_011;
            rec4 = [path_rec{2}, 1];
        else
            weight(4, i) = w_101_011;
            rec4 = [path_rec{6}, 1];
        end

        if w_010_100 > w_110_100
            weight(5, i) = w_010_100;
            rec5 = [path_rec{3}, 0];
        else
            weight(5, i) = w_110_100;
            rec5 = [path_rec{7}, 0];
        end

        if w_010_101 > w_110_101
            weight(6, i) = w_010_101;
            rec6 = [path_rec{3}, 1];
        else
            weight(6, i) = w_110_101;
            rec6 = [path_rec{7}, 1];
        end

        if w_011_110 > w_111_110
            weight(7, i) = w_011_110;
            rec7 = [path_rec{4}, 0];
        else
            weight(7, i) = w_111_110;
            rec7 = [path_rec{8}, 0];
        end

        if w_011_111 > w_111_111
            weight(8, i) = w_011_111;
            rec8 = [path_rec{4}, 1];
        else
            weight(8, i) = w_111_111;
            rec8 = [path_rec{8}, 1];
        end

        path_rec{1} = rec1;
        path_rec{2} = rec2;
        path_rec{3} = rec3;
        path_rec{4} = rec4;
        path_rec{5} = rec5;
        path_rec{6} = rec6;
        path_rec{7} = rec7;
        path_rec{8} = rec8;
    end
end

[~, idx] = max(weight(:, nn));
out = path_rec{idx};
    
end