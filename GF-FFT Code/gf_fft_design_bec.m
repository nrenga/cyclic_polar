function [ f ] = gf_fft_design_bec(n,e,d,q)
% Design FFT code of length N=2^n for BEC(e) and target block error rate d
% Over GF(q) where q is prime
% Generate virtual channel erasure probabilities
% f = input a priori probs matrix in input order;
%     each column signifies a particular bit; each row signifies a GF value
E = e;
for i = 1:n
    % Interleave updates to keep in FFT decoding order
    E = reshape([1-(1-E).*(1-E); E.*E],1,[]);
end
N = 2^n;
bin_fft_order = fliplr(dec2bin(0:(N-1),log2(N)));
fft_order = bin2dec(bin_fft_order) + 1;
E = E(fft_order);
% Sort into increasing order and compute cumulative sum
[SE,order] = sort(E);
CSE = cumsum(SE);
% Find good indices
I = sum(double(CSE<d));
f = zeros(q,2^n);
f(:,order(1:I)) = 1/q;    % For good bits, all 'q' values are equi-probable
f(1,f(1,:)==0) = 1;       % Set frozen bits to 0s
end