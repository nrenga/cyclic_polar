function [res] = gfpow(a,b,field)
% Computes a^b in GF(p^m) represented as field
res = 0;   % 0 is 1 in exponential format
for times = 1:b
    res = gfmul(res,a,field);
end
return
end