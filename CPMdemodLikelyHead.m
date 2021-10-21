function out = CPMdemodLikelyHead(CPM_recv, oversamp, nn)
%CPMdemodLikely - Description
%
% Syntax: out = CPMdemodLikely(CPM_recv)
%
% Long description
    out = zeros(1, nn);
    % weight2 = zeros(2, 1);
    weight4 = zeros(4, 1);
    weight4End = zeros(4, 1);
    weight8 = zeros(8, 1);
    sample = 4;
    code = [[0,0,0]; [0,1,0]; [1,0,0]; [1,1,0]];
    % code = [[0,0,0]; [0,0,1]; [0,1,0]; [0,1,1]; [1,0,0]; [1,0,1]; [1,1,0]; [1,1,1]];

    phi = zeros(4, 1);
    % phi2 = zeros(2, 1);
    phi_00 = exp(-1i*((1:sample)-1)*pi/2/sample);
    phi_01 = ones(1, sample);
    phi_10 = ones(1, sample);
    phi_11 = exp(1i*((1:sample)-1)*pi/2/sample);
    % phi_100 = [ones(1, sample), exp(-1i*((1:sample)-1)*pi/2/sample)];
    % phi_101 = ones(1, sample*2);
    % phi_110 = [exp(1i*((1:sample)-1)*pi/2/sample), exp(1i*pi/2)*ones(1, sample)];
    % phi_111 = exp(1i*((1:sample*2)-1)*pi/sample/2);

    %  pi/4 100 101 110 111
    phi_00_s = [ones(1, sample), exp(-1i*((1:sample)-1)*pi/2/sample)];
    phi_01_s = [ones(1, sample), ones(1, sample)];
    phi_10_s = [exp(1i*((1:sample)-1)*pi/2/sample), exp(1i*pi/2)*ones(1, sample)];
    phi_11_s = [exp(1i*((1:sample)-1)*pi/2/sample), exp(1i*pi/2)*exp(1i*((1:sample)-1)*pi/2/sample)];

    for i = 2:nn
        CPM = CPM_recv((i-1)*sample+1:i*sample);
        if i == 2
            CPM = CPM_recv((i-2)*sample+1:i*sample);
            weight4(1) = real(sum(CPM .* conj(phi_00_s)));
            weight4(2) = real(sum(CPM .* conj(phi_01_s)));
            weight4(3) = real(sum(CPM .* conj(phi_10_s)));
            weight4(4) = real(sum(CPM .* conj(phi_11_s)));
            phi = [exp(-1i*pi/2); 1; exp(1i*pi/2); exp(1i*pi)];
        elseif i == nn
            weight4End(1) = weight4(1) + real(sum(CPM .* conj(phi(1) .* phi_00)));
            weight4End(2) = weight4(2) + real(sum(CPM .* conj(phi(2) .* phi_10)));
            weight4End(3) = weight4(3) + real(sum(CPM .* conj(phi(3) .* phi_00)));
            weight4End(4) = weight4(4) + real(sum(CPM .* conj(phi(4) .* phi_10)));
            [~, idx] = max(weight4End);
            out(i-2:end) = code(idx);
        else
            weight8(1) = weight4(1) + real(sum(CPM .* conj(phi(1) .* phi_00)));
            weight8(2) = weight4(1) + real(sum(CPM .* conj(phi(1) .* phi_01)));
            weight8(3) = weight4(2) + real(sum(CPM .* conj(phi(2) .* phi_10)));
            weight8(4) = weight4(2) + real(sum(CPM .* conj(phi(2) .* phi_11)));
            weight8(5) = weight4(3) + real(sum(CPM .* conj(phi(3) .* phi_00)));
            weight8(6) = weight4(3) + real(sum(CPM .* conj(phi(3) .* phi_01)));
            weight8(7) = weight4(4) + real(sum(CPM .* conj(phi(4) .* phi_10)));
            weight8(8) = weight4(4) + real(sum(CPM .* conj(phi(4) .* phi_11)));  
            phi_tmp = [phi(1)*exp(-1i*pi/2); phi(1); phi(2); phi(2)*exp(1i*pi/2); phi(3)*exp(-1i*pi/2); phi(3); phi(4); phi(4)*exp(1i*pi/2)];
            [~, idx] = max(weight8);
            % if i == nn
            %     out(i-2:end) = code(idx,:);
            % else
                if idx <= 4
                    out(i-2) = 0;
                    weight4 = weight8(1:4);
                    % 0011
                    phi = phi_tmp(1:4);
                else
                    out(i-2) = 1;
                    weight4 = weight8(5:8);
                    % 0011
                    phi = phi_tmp(5:8);
                end
            % end
        end
    end
    out = out(1:end-2);
end