clc
close all
clear all
%N = input('Give the blocklength: N = ');
Nmax = 5000;
prime_len = primes(Nmax);
%prime_len = 2.^[1:16];
max_rates = zeros(length(prime_len),2);
for iter = 1:length(prime_len)
    N = prime_len(iter);
    e = 0.5;
    d = 0.1;
    q = 7;   % Doesn't matter for DE
    lengths = factor(N);
    fprintf('Blocklength = %d\n',prod(lengths));
    if (all(lengths == lengths(1)))
        I = lengths(1,:);
    else
        I = perms(lengths);
    end
    k = zeros(size(I,1),1);
    for i = 1:size(I,1)
        f = gf_fft_design_qec(I(i,:),e,d,q);
        k(i) = nnz(f == 1/q);
    end
    unique_rates = unique(k);
    unique_rates(:,2) = unique_rates(:,1)/N;
    fprintf('\n\tk\tRate\n');
    disp(unique_rates);
    max_rates(iter,1) = N;
    max_rates(iter,2) = max(unique_rates(:,2));
end
figure;
plot(max_rates(:,1),max_rates(:,2));
xlabel('Blocklength N','FontSize',16);
ylabel('Rate = k/N','FontSize',16);
title('Max. Rates of GF-FFT Code - Prime Lengths; e = 0.5; Blockerr = 0.1','FontSize',16);