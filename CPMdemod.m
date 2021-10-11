function out = CPMdemod(CPM_recv)
%CPMdemod - Description
%
% Syntax: out = CPMdemod(CPM_recv)
%
% Long description
nn = 50;
path_rec = cell(4, 1);
weight = -10*ones(4, nn);
phi = zeros(4, 1);

down_sample = 4;

pi_4_p = exp(1i * ((1:down_sample)-1) * (pi/4/down_sample));
pi_4_n = exp(-1i * ((1:down_sample)-1) * (pi/4/down_sample));
pi_2_p = exp(1i * ((1:down_sample)-1) * (pi/2/down_sample));
pi_2_n = exp(-1i * ((1:down_sample)-1) * (pi/2/down_sample));

for i = 1:nn
    if i == 1
        weight(1, i) = real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(pi_4_p)));
        weight(2, i) = real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(pi_4_n)));
        path_rec{1} = [0];
        path_rec{2} = [1];
        phi(1) = exp(-1i*pi/4);
        phi(2) = exp(1i*pi/4);
    else
        rec1 = [];
        rec2 = [];
        rec3 = [];
        rec4 = [];
        phi1 = [];
        phi2 = [];
        phi3 = [];
        phi4 = [];
        %% 00
        % 0
        % w1_0 = weight(1, i-1) + real(CPM_recv(i)*conj(phi(1)*exp(-1i*pi/2)));
        w1_0 = weight(1, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(1)*pi_2_n)));
        % 1
        w1_1 = weight(1, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(1))));

        %% 10
        % 0
        % w3_0 = weight(3, i-1) + real(CPM_recv(i)*conj(phi(3)*exp(-1i*pi/2)));
        w3_0 = weight(3, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(3)*pi_2_n)));
        if w3_0 > w1_0
            weight(1, i) = w3_0;
            rec1 = [path_rec{3}, 0];
            phi1 = phi(3)*exp(-1i*pi/2);
        else
            weight(1, i) = w1_0;
            rec1 = [path_rec{1}, 0];
            phi1 = phi(1)*exp(-1i*pi/2);
        end
        % 1
        w3_1 = weight(3, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(3))));
        if w3_1 > w1_1
            weight(2, i) = w3_1;
            rec2 = [path_rec{3}, 1];
            phi2 = phi(3);
        else
            weight(2, i) = w1_1;
            rec2 = [path_rec{1}, 1];
            phi2 = phi(1);
        end

        %% 01
        % 0
        w2_0 = weight(2, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(2))));
        % 1
        % w2_1 = weight(2, i-1) + real(CPM_recv(i)*conj(phi(2)*exp(1i*pi/2)));  
        w2_1 = weight(2, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(2)*pi_2_p)));    
        
        %% 11
        % 0
        w4_0 = weight(4, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(4))));
        if w4_0 > w2_0
            weight(3, i) = w4_0;
            rec3 = [path_rec{4}, 0];
            phi3 = phi(4);
        else
            weight(3, i) = w2_0;
            rec3 = [path_rec{2}, 0];
            phi3 = phi(2);
        end
        % 1
        % w4_1 = weight(4, i-1) + real(CPM_recv(i)*conj(phi(4)*exp(1i*pi/2)));
        w4_1 = weight(4, i-1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(4)*pi_2_p)));
        if w4_1 > w2_1
            weight(4, i) = w4_1;
            rec4 = [path_rec{4}, 1];
            phi4 = phi(4)*exp(1i*pi/2);
        else
            weight(4, i) = w2_1;
            rec4 = [path_rec{2}, 1];
            phi4 = phi(2)*exp(1i*pi/2);
        end
        path_rec{1} = rec1;
        path_rec{2} = rec2;
        path_rec{3} = rec3;
        path_rec{4} = rec4;
        phi = [phi1; phi2; phi3; phi4];
    end
end

[~, idx] = max(weight(:, nn));
out = path_rec{idx};
    
end