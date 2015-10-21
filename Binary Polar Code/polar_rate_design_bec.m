function f = polar_rate_design_bec(n,k)
% Design polar code of length N=2^n for BEC(e) and target block error rate d
% Generate virtual channel erasure probabilities
E = 0.5;   % Design Erasure Rate
for i=1:n
    % Interleave updates to keep in polar decoding order
    E = reshape([1-(1-E).*(1-E); E.*E],1,[]);
end
% Sort into increasing order and compute cumulative sum
[SE,order] = sort(E);
f = zeros(1,length(E));
f(order(1:k)) = 1/2;    % For good bits, "1" and "0" are equi-probable. So P(1) = 1/2
end