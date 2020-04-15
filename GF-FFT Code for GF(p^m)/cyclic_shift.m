function [x] = cyclic_shift(u, s)

len = length(u);

if (size(u,1) == 1)
    x = [u((len-s+1):len) u(1:(len-s))];
else
    x = [u((len-s+1):len) ; u(1:(len-s))];
end

end