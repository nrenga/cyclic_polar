function [ f ] = cp_rate_design_qec(n,k,q)
% Design FFT code of length N=2^n for BEC(e) and target block error rate d
% Over GF(q) where q is prime
% Generate virtual channel erasure probabilities
% f = input a priori probs matrix in input order;
%     each column signifies a particular bit; each row signifies a GF value
E = 0.5;    % Design Erasure Rate
for i = 1:n
    % Interleave updates to keep in FFT decoding order;
    E = reshape([1-(1-E).*(1-E); E.*E],1,[]);
end
N = 2^n;
bin_fft_order = fliplr(dec2bin(0:(N-1),log2(N)));
fft_order = bin2dec(bin_fft_order)+1;
E = E(fft_order);
% Sort into increasing order and compute cumulative sum
[SE,order] = sort(E);
% Find good indices
f = zeros(q,N);
f(:,fft_order(order(1:k))) = 1/q;    % For good bits, all 'q' values are equi-probable
f(1,f(1,:)==0) = 1;       % Set frozen bits to 0s
end