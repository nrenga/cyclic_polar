% Code to test Forney for FFT code
clc
close all
clear all

N = 18;
p = 19;
m = 1;
q = p^m;
% verbose = 0 -> no intermediate outputs are displayed
% verbose = 1 -> sufficient intermediate outputs are displayed
% verbose = 2 -> most intermediate outputs are displayed
verbose = 0;  

% a = 6; b = 3;
% inds = 1 + reshape(reshape(0:(N-1),a,b)',1,[]);
% ins2 = sub2ind([N,N],1:N,inds);
% P = zeros(N);
% P(ins2) = 1;
% Pretro = [1 zeros(1,N-1); zeros(N-1,1) fliplr(eye(N-1))];
% % P2 = P*diag(2*(inds-1));
% % P2(P2==0) = -Inf; P2(1,1) = 0;
% Pf = P; Pf(Pf==0) = -Inf; Pf(Pf==1) = 0;

prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);
    
if (mod(p^m-1,N) == 0)
    lengths = factor(N);
    e = 0.5;
    d = 0.1;
    [f, I, O] = gf_fft_design_qec(lengths,e,d,q);
    k = nnz(f == 1/q);
    
    % No of errors/erasures that can be corrected
    %t = floor((N-k)/2);
    
    u = (-Inf)*ones(1,N);  % In exponential format;
    u(f == 1/q) = randi(q,1,k)-2;
    u(u == -1) = -Inf;
%     u = [ -Inf  -Inf  -Inf  -Inf  -Inf  -Inf  -Inf  -Inf    14  -Inf  -Inf  -Inf  -Inf     6     5];
%     u = [ -Inf  -Inf  -Inf  -Inf  -Inf  -Inf  -Inf  -Inf     4  -Inf  -Inf     8  -Inf     2    11];
    
    fprintf('\nk = %d, N = %d, q = %d^%d = %d, Rate = %1.4f.\n\nEncoding now...\n',k,N,p,m,q,k/N);
    tic
    [v] = gf_ifft(u,p,m,[],N);
    toc
    %[v, Winv] = naive_gf_ifft(u,p,m);
    
    [v2, W] = naive_gf_fft(u,p,m);
%     v = [4     8     0     6     6     5     3     0     0    10     8    13  -Inf    13     7]';
%     v = [ 6    13     7    13     3     9     7     4    13     1     5     5     3  -Inf     9]';
    %fprintf('\nCodeword:');
    %disp(v);
    
    r = v;
    
    % QSCE Channel
    perror = 0;
    perasure = 0.55;
    pcorrect = 1 - perror - perasure;
    symbols = [1/q -Inf 0:(q-2)];
    for i = 1:N
        pvec = [perasure (perror/(q-1))*ones(1,q)];
        pvec(symbols == v(i)) = pcorrect;
        plen = length(pvec);     % plen = q+1
        pvec = reshape(pvec,plen,1);   % pvec = pvec'
        index = sum(repmat(rand(1),plen,1) > cumsum(pvec),1) + 1;
        r(i) = symbols(index);
    end
%     r = [4     8     0     6     6     5     3     0     0    10     8     8  -Inf    13     8]';
%     r = [6    13     7    13     8     9     4     4    13     1     5     5     3     4     9]';
    
    erased_indices = find(r == 1/q);
    fprintf('\nNumber of erasures in received word = %d\n',nnz(erased_indices));
    
    unerased_indices = 1:N;
    unerased_indices(erased_indices) = [];
    
    error_indices = unerased_indices(r(unerased_indices) ~= v(unerased_indices));
    fprintf('\nNumber of errors in received word = %d\n',nnz(error_indices));
    
%     fprintf('\nFor the largest block in the code with the maximum number of roots, in the worst case we have:\n');
%     fprintf('\ndmin = %d; 2*no_of_errors + no_of_erasures + 1 = %d\n',max(lengths)+1,2*nnz(error_indices)+nnz(erased_indices)+1);
    
%     if (max(lengths)+1 < 2*nnz(error_indices)+nnz(erased_indices)+1)
%         fprintf('\nDECODER FAILURE POSSIBLE!\n');
%     end
   
    % Now decode
    tic
    [u_est,r_est] = gf_fft_decode_qsce(r,f,p,m,[],lengths,N,verbose);
    toc
    u_est = u_est';
    %fprintf('\nTx:');
    %disp(v);
    %fprintf('Rx:');
    %disp(r);
    %fprintf('Est:');
    %disp(r_est');
        
    fprintf('\nNumber of decoding errors = %d\n\n',nnz(r_est ~= v));
else
    fprintf('\nN = %d does not divide p^m-1 = %d. Check.\n',N,p^m-1);
end

% % Permutation symmetry of DFT
% perm_indices2 = mod(2*[0:8],9);
% mat_inds2 = sub2ind(size(W),1:9,perm_indices2+1);
% P2 = zeros(9); P2(mat_inds2) = 1;
% perm_indices = mod(5*[0:8]+2,9);
% mat_inds = sub2ind(size(W),1:9,perm_indices+1);
% P1 = zeros(9); P1(mat_inds) = 1;
% D  = diag(mod(-(2)*(2*2)*[0:8],18));  % omega = alpha^2
% D(D==0) = -Inf; D(1,1) = 0;
% D  = diag(mod(-(2)*(2*2)*[0:8],18));
% DP2 = D*P2;
% DP2(DP2==0) = -Inf; DP2(1,1) = 0;
% 
% % Perfect-shuffle permutation induced linear transform matrix M
% Pretro2 = Pretro; Pretro2(Pretro2==0) = -Inf; Pretro2(Pretro2==1) = 0;
% n = -Inf;       % -1 is 0 in exponential format
% for i = 1:N
% n = gfadd(n,0,field);  % Compute N in GF(p^m)
% end
% Ninv = gfdiv(0,n,field);
% A = Ninv*eye(N); A(A==0) = -Inf;
% M = gf_mat_multiply(A,gf_mat_multiply(Pretro2,gf_mat_multiply(W',P*W,field),field),field);
