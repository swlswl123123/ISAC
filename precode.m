function out = precode(bits)
%precode - Description
%
% Syntax: out = precode(bits)
%
% Long description
    out = zeros(size(bits));
    for i = 1:length(bits)
        if rem(i, 2)
            if i == 1
                out(i) = xor(bits(i), 0);
            else
                out(i) = xor(bits(i), bits(i-1)); 
            end
        else
            out(i) = ~xor(bits(i), bits(i-1));
        end
    end
end