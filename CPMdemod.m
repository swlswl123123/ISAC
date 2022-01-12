function out = CPMdemod(CPM_recv, down_sample, nn)
%CPMdemod - Description
%
% Syntax: out = CPMdemod(CPM_recv)
%
% Long description

path_rec = cell(4, 1);
weight = zeros(4, 1);
phi = zeros(4, 1);

pi_2_p = exp(1i * (1:down_sample) * (pi/2/down_sample));
pi_2_n = exp(-1i * (1:down_sample) * (pi/2/down_sample));

for i = 1:nn
    if i == 1
        weight(1) = real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(pi_2_n)));
        weight(2) = real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample)));
        phi(1) = exp(-1i*pi/2);
        phi(2) = exp(0);
    elseif i == 2
        w_00_00 = weight(1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(1)*pi_2_n)));
        phi_00_00 = phi(1) * exp(-1i*pi/2);
        w_00_01 = weight(1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(1))));
        phi_00_01 = phi(1);
        w_01_10 = weight(2) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(2))));
        phi_01_10 = phi(2);
        w_01_11 = weight(2) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(2)*pi_2_p)));
        phi_01_11 = phi(2) * exp(1i*pi/2);

        weight(1) = w_00_00;
        weight(2) = w_00_01;
        weight(3) = w_01_10;
        weight(4) = w_01_11;

        phi(1) = phi_00_00;
        phi(2) = phi_00_01;
        phi(3) = phi_01_10;
        phi(4) = phi_01_11;

        path_rec{1} = [0, 0];
        path_rec{2} = [0, 1];
        path_rec{3} = [1, 0];
        path_rec{4} = [1, 1];

    else
        w_00 = [];
        w_01 = [];
        w_10 = [];
        w_11 = [];

        phi_00 = [];
        phi_01 = [];
        phi_10 = [];
        phi_11 = [];

        rec_00 = [];
        rec_01 = [];
        rec_10 = [];
        rec_11 = [];

        % 00
        w_00_00 = weight(1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(1)*pi_2_n)));
        w_10_00 = weight(3) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(3)*pi_2_n)));
        if w_00_00 > w_10_00
            w_00 = w_00_00;
            rec_00 = [path_rec{1}, 0];
            phi_00 = phi(1) * exp(-1i*pi/2);
        else
            w_00 = w_10_00;
            rec_00 = [path_rec{3}, 0];
            phi_00 = phi(3) * exp(-1i*pi/2);   
        end

        % 01
        w_00_01 = weight(1) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(1))));
        w_10_01 = weight(3) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(3))));
        if w_00_01 > w_10_01
            w_01 = w_00_01;
            rec_01 = [path_rec{1}, 1];
            phi_01 = phi(1);
        else
            w_01 = w_10_01;
            rec_01 = [path_rec{3}, 1];
            phi_01 = phi(3);
        end

        % 10
        w_01_10 = weight(2) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(2))));
        w_11_10 = weight(4) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(4))));
        if w_01_10 > w_11_10
            w_10 = w_01_10;
            rec_10 = [path_rec{2}, 0];
            phi_10 = phi(2);
        else
            w_10 = w_11_10;
            rec_10 = [path_rec{4}, 0];
            phi_10 = phi(4);
        end

        % 11
        w_01_11 = weight(2) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(2)*pi_2_p)));
        w_11_11 = weight(4) + real(sum(CPM_recv((i-1)*down_sample+1:i*down_sample).*conj(phi(4)*pi_2_p)));    
        if w_01_11 > w_11_11
            w_11 = w_01_11;
            rec_11 = [path_rec{2}, 1];
            phi_11 = phi(2) * exp(1i*pi/2);
        else
            w_11 = w_11_11;
            rec_11 = [path_rec{4}, 1];
            phi_11 = phi(4) * exp(1i*pi/2);
        end

        weight = [w_00; w_01; w_10; w_11];

        phi = [phi_00; phi_01; phi_10; phi_11];

        path_rec{1} = rec_00;
        path_rec{2} = rec_01;
        path_rec{3} = rec_10;
        path_rec{4} = rec_11;
    end
end

[~, idx] = max(weight);
out = path_rec{idx};
    
end