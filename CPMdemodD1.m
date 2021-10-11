function out = CPMdemodD1(CPM_recv)
%CPMdemod - Description
%
% Syntax: out = CPMdemod(CPM_recv)
%
% Long description
nn = 50;
path_rec = cell(4, 1);
weight = -10*ones(4, nn);

down_sample = 4;

pi_00_00 = exp(-1i * pi/2) * ones(1, down_sample);
pi_00_01 = exp(1i * (((1:down_sample)-1) * (pi/2/down_sample) - pi/2));
pi_01_10 = ones(1, down_sample);
pi_01_11 = exp(1i * ((1:down_sample)-1) * (pi/2/down_sample));
pi_10_00 = exp(-1i * ((1:down_sample)-1) * (pi/2/down_sample));
pi_10_01 = ones(1, down_sample);
pi_11_10 = exp(1i * (-((1:down_sample)-1) * (pi/2/down_sample) + pi/2));
pi_11_11 = exp(1i * pi/2) * ones(1, down_sample);

for i = 2:nn
    CPM_dif = CPM_recv((i-1)*down_sample+1:i*down_sample) .* conj(CPM_recv((i-2)*down_sample+1:(i-1)*down_sample));
    if i == 2
        pi_00_00_t = exp(1i * (-((1:down_sample)-1) * (pi/4/down_sample) - pi/4));
        pi_00_01_t = exp(1i * (((1:down_sample)-1) * (pi/4/down_sample) - pi/4));
        pi_11_10_t = exp(1i * (-((1:down_sample)-1) * (pi/4/down_sample) + pi/4));
        pi_11_11_t = exp(1i * (((1:down_sample)-1) * (pi/4/down_sample) + pi/4));
        weight(1, i) = real(sum(CPM_dif.*conj(pi_00_00_t)));
        weight(2, i) = real(sum(CPM_dif.*conj(pi_00_01_t)));
        weight(3, i) = real(sum(CPM_dif.*conj(pi_11_10_t)));
        weight(4, i) = real(sum(CPM_dif.*conj(pi_11_11_t)));
        path_rec{1} = [0, 0];
        path_rec{2} = [0, 1];
        path_rec{3} = [1, 0];
        path_rec{4} = [1, 1];
    else
        rec1 = [];
        rec2 = [];
        rec3 = [];
        rec4 = [];
        %% 00
        % 0
        w1_0 = weight(1, i-1) + real(sum(CPM_dif.*conj(pi_00_00)));
        % 1
        w1_1 = weight(1, i-1) + real(sum(CPM_dif.*conj(pi_00_01)));

        %% 10
        % 0
        w3_0 = weight(3, i-1) + real(sum(CPM_dif.*conj(pi_10_00)));
        if w3_0 > w1_0
            weight(1, i) = w3_0;
            rec1 = [path_rec{3}, 0];
        else
            weight(1, i) = w1_0;
            rec1 = [path_rec{1}, 0];
        end
        % 1
        w3_1 = weight(3, i-1) + real(sum(CPM_dif.*conj(pi_10_01)));
        if w3_1 > w1_1
            weight(2, i) = w3_1;
            rec2 = [path_rec{3}, 1];
        else
            weight(2, i) = w1_1;
            rec2 = [path_rec{1}, 1];
        end

        %% 01
        % 0
        w2_0 = weight(2, i-1) + real(sum(CPM_dif.*conj(pi_01_10)));
        % 1
        w2_1 = weight(2, i-1) + real(sum(CPM_dif.*conj(pi_01_11)));    
        
        %% 11
        % 0
        w4_0 = weight(4, i-1) + real(sum(CPM_dif.*conj(pi_11_10)));
        if w4_0 > w2_0
            weight(3, i) = w4_0;
            rec3 = [path_rec{4}, 0];
        else
            weight(3, i) = w2_0;
            rec3 = [path_rec{2}, 0];
        end
        % 1
        w4_1 = weight(4, i-1) + real(sum(CPM_dif.*conj(pi_11_11)));
        if w4_1 > w2_1
            weight(4, i) = w4_1;
            rec4 = [path_rec{4}, 1];
        else
            weight(4, i) = w2_1;
            rec4 = [path_rec{2}, 1];
        end
        path_rec{1} = rec1;
        path_rec{2} = rec2;
        path_rec{3} = rec3;
        path_rec{4} = rec4;
    end
end

[~, idx] = max(weight(:, nn));
out = path_rec{idx};
    
end