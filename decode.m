function out = decode(code)
%decode - Description
%
% Syntax: out = decode(code)
%
% Long description
    out = zeros(size(code));
    for i = 1:length(code)
        if rem(i, 2)
            if i == 1
                out(i) = xor(0, code(i));
            else
                out(i) = xor(code(i), out(i-1));
            end
        else
            out(i) = ~xor(code(i), out(i-1));
        end
    end
end