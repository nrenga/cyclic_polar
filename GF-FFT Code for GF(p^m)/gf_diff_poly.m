function [px_d] = gf_diff_poly(px,field)
% Differentiate the polynomial px in GF(p^m)

px_d = px(2:end);
for k = 1:length(px_d)
    multiplier = gf_find_n(k,field);
    px_d(k) = gfmul(px_d(k),multiplier,field);
end
return
end