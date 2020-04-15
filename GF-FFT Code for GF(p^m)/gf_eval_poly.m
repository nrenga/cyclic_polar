function [val] = gf_eval_poly(px,point,field)
% Evaluate the polynomial px at the value 'point' in GF(p^m)
% 'point' is in exponential format
% px is a vector containing coefficients of p(x) in ascending order of
% exponents of x i.e. p(x)=px(1)+px(2)*x+px(3)*x^2+...
% The components of px are in exponential format too

val = -1;
for k = 0:(length(px)-1)
    multiplier = gfpow(point,k,field);
    val = gfadd(val,gfmul(px(k+1),multiplier,field),field);
end
return
end