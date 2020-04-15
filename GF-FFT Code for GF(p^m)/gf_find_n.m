function [res] = gf_find_n(n,field)
% Find the value of an integer 'n' in GF(p^m)
% n = 1+1+1+... (n times)
% Returned value is in exponential format

res = -1;       % -1 is 0 in exponential format
for times = 1:n
    res = gfadd(res,0,field);  % Compute n in GF(p^m)
end
return
end