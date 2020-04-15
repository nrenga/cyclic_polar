% Code to test Forney for FFT code over a channel with errors & erasures
clc
close all
clear all

N = 15;
p = 2;
m = 4;
q = p^m;
prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);
lengths = factor(N);
e = 0.5;
d = 0.1;
f = gf_fft_design_qec(lengths,e,d,q);
% k = nnz(f == 1/q);
k = 9;
verbose = 2;

% No of errors/erasures that can be corrected
t = floor((N-k)/2);

b = 0;      % b <= t+1
if (b <= k)
    %u = [randi(q,1,b)-2 (-Inf)*ones(1,N-k) randi(q,1,k-b)-2];  % In exponential format; alpha^9 = 3
    u = [5 7 1 -Inf 4 10 9 7 1 3 13 1 12 5 1]; 
    
    fprintf('\nk = %d, N = %d, q = %d^%d = %d, Rate = %1.4f.\n\nEncoding now...\n',k,N,p,m,q,k/N);
    v = gf_ifft(u,p,m,[],N);
    %v = [ 10     1     3    14    -1     2     3    13     1    14     7     5    12     2     3];
    
    %fprintf('\nCodeword:');
    %disp(v);
    
    %r = v;
    
    r = [2 1 3 5 8 10 7 5 13 6 3 6 10 0 7]';
        
    % QSCE Channel
%     perror = 0;
%     perasure = 0.1;
%     pcorrect = 1 - perror - perasure;
%     symbols = [1/q -Inf 0:(q-2)];
%     for i = 1:N
%         pvec = [perasure (perror/(q-1))*ones(1,q)];
%         pvec(symbols == v(i)) = pcorrect;
%         plen = length(pvec);     % plen = q+1
%         pvec = reshape(pvec,plen,1);   % pvec = pvec'
%         index = sum(repmat(rand(1),plen,1) > cumsum(pvec),1) + 1;
%         r(i) = symbols(index);
%     end
    
    % Introduce errors
%     err = 0.1;
%     error_indices = find(rand(1,N)<err);
%     for j = 1:length(error_indices)
%         set = [-Inf 0:(q-2)];
%         set(set==r(error_indices(j))) = [];
%         r(error_indices(j)) = set(randi(q-1));
%     end    
    
    % Introduce erasures
    %err = t/N;
%     err = 0.3;
%     erased_indices = find(rand(1,N)<err);
%     %erased_indices = [3 4 8 10 14];
%     r(erased_indices) = 1/q;  % Symbols erased

    erased_indices = find(r == 1/q);
    fprintf('\nNumber of erasures in received word = %d\n',nnz(erased_indices));
    
    unerased_indices = 1:N;
    unerased_indices(erased_indices) = [];
    
    error_indices = unerased_indices(r(unerased_indices) ~= v(unerased_indices));
    fprintf('\nNumber of errors in received word = %d\n',nnz(error_indices));
    
    fprintf('\ndmin = %d; 2*no_of_errors + no_of_erasures + 1 = %d\n',N-k+1,2*nnz(error_indices)+nnz(erased_indices)+1);
    
    if (N-k+1 < 2*nnz(error_indices)+nnz(erased_indices)+1)
        fprintf('\nDECODER FAILURE!\n');
    end

    %inp = [(1/q)*ones(1,b) (-Inf)*ones(1,N-k) (1/q)*ones(1,k-b)];   % Known information at input side
    inp = [u(1:6) (1/q)*ones(1,9)];
    
    % Now decode
    [u_est,r_est] = bm_forney_fft(r,inp,p,m,verbose);
    u_est = u_est';
    %fprintf('\nTx:');
    %disp(v);
    %fprintf('Rx:');
    %disp(r);
    %fprintf('Est:');
    %disp(r_est');
        
    fprintf('\nNumber of decoding errors = %d\n\n',nnz(r_est ~= v));
else
    fprintf('\nFor the code, k=%d and so, b should be at most k. But, b=%d',k,b);
end