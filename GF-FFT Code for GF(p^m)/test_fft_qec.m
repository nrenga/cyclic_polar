% Code to test Forney for FFT code
clc
close all
clear all

N = 30;
p = 31;
m = 1;
q = p^m;

if (mod(p^m-1,N) == 0)
    lengths = factor(N);
    e = 0.5;
    d = 0.1;
    f = gf_fft_design_qec(lengths,e,d,q);
    k = nnz(f == 1/q);
    
    % No of errors/erasures that can be corrected
    %t = floor((N-k)/2);
    
    u = (-Inf)*ones(1,N);  % In exponential format;
    u(f == 1/q) = randi(q,1,k)-2;
    u(u == -1) = -Inf;
    %u = [ -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 27 -Inf -Inf 23 -Inf 14 12];
    
    fprintf('\nk = %d, N = %d, q = %d^%d = %d, Rate = %1.4f.\n\nEncoding now...\n',k,N,p,m,q,k/N);
    v = gf_ifft(u,p,m,[],N);
    %fprintf('\nCodeword:');
    %disp(v);
    
    % Receiver
    r = v;
    err = 0;
    %erased_indices = [1 2 8 11];
    %r(erased_indices) = 1/q;
    r(rand(1,N)<err) = 1/q;  % Symbols erased
    fprintf('\nNumber of erasures in received word = %d\n',nnz(r == 1/q));
    
    fprintf('\nStarting to decode...\n');
    [u_est,r_est] = gf_fft_decode_qec(r,f,p,m,[],lengths,N);
    u_est = u_est';
    %fprintf('Actual Input:');
    %disp(u);
    %fprintf('Decoded Input:');
    %disp(inp_est);
    fprintf('\nNumber of decoding errors = %d\n\n',nnz(v~=r_est));
else
    fprintf('\nN = %d does not divide p^m-1 = %d. Check.\n',N,p^m-1);
end