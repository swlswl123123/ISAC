function out = pn_gen(data)
%pn_gen - Description
%
% Syntax: out = pn_gen(data)
%
% Long description
    data = 2*data - 1;
    N = length(data);
    data = reshape(data, [2, N / 2]);
    I = repmat(data(2, :), [2 1]);
    I = reshape(-I, [1, N]);
    for i = 1:N
        if i > 1 && I(i-1) * I(i) < 0
            I(i) = 0;
        end
    end
    I = I(1:end-1);
    Q = repmat(data(1, :), [2 1]);
    Q = reshape(Q, [1 N]);
    for i = 1:N
        if i > 1 && Q(i-1) * Q(i) < 0
            Q(i) = 0;
        end
    end
    Q = Q(2:end);
    out = I + 1i*Q;
end