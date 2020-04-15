% Code to test Forney for FFT code
clc
close all
clear all

N = 15;
p = 31;
m = 1;
q = p^m;
prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);
lengths = factor(N);
e = 0.5;
d = 0.1;
f = gf_fft_design_qec(lengths,e,d,q);
k = nnz(f == 1/q);

% No of errors/erasures that can be corrected
t = floor((N-k)/2);

b = 0;      % b <= t+1
if (b <= k)
    u = [randi(q,1,b)-2 (-1)*ones(1,N-k) randi(q,1,k-b)-2];  % In exponential format; alpha^9 = 3
    
    fprintf('\nk = %d, N = %d, q = %d^%d = %d, Rate = %1.4f.\n\nEncoding now...\n',k,N,p,m,q,k/N);
    v = gf_fft(u,p,m,[],N);
    v(v == -Inf) = -1;
    %fprintf('\nCodeword:');
    %disp(v);
    
    % Receiver
    r = v;
    err = t/N;
    r(rand(1,N)<err) = 1/q;  % Symbols erased
    fprintf('\nNumber of erasures in received word = %d\n',nnz(r == 1/q));
    
    inp = [(1/q)*ones(1,b) (-1)*ones(1,N-k) (1/q)*ones(1,k-b)];   % Known information at input side
    
    [u_est,r_est] = forney_fft(r,inp,p,m);
    u_est = u_est';
    %fprintf('Actual Input:');
    %disp(u);
    %fprintf('Decoded Input:');
    %disp(inp_est);
    
    fprintf('Number of decoding errors = %d\n\n',nnz(v-r_est'));
else
    fprintf('\nFor the code, k=%d and so, b should be at most k. But, b=%d',k,b);
end