clc
close all
clear all
t1 = clock;
% Setup code parameters
N = 255;
e = 0.5;   % Design channel erasure rate
d = 0.1;   % Design block erasure rate
p = 2;
m = 8;
q = p^m;

if (mod(q-1,N) == 0)
    lengths = factor(N);
    f = gf_fft_design_qec(lengths,e,d,q);  % Design is same with and without twiddle
    k = nnz(f == 1/q);
    rate = k/N;
    
    % Find best frozen bits and compute rate
    % f = input a priori probs matrix in input order;
    %     each column signifies a particular bit; each row signifies a GF value
    M = 100;
    blockerr = [];
    for e = 0.1:0.1:1
        % Set frozen bits, add random data, and encode
        u = (-1)*ones(1,N);  % In exponential format;
        u(f == 1/q) = randi(q,1,k)-2;
        fprintf('\nk = %d, N = %d, q = %d^%d = %d, Rate = %1.4f.\n\nEncoding now...\n',k,N,p,m,q,rate);
        x = gf_fft(u,p,m,[],N);
        
        % Run a few sims
        biterr = zeros(1,M);
        for i = 1:M
            clc
            fprintf('\nErasure Rate: %f\n',e);
            fprintf('\nSimulation %d out of %d\n',i,M);
            disp(blockerr);
            
            % Transmit
            recv = x;
            recv(rand(1,N)<e) = 1/q;  % Random erasures
            fprintf('\nNumber of erasures in received word = %d\n',nnz(recv == 1/q));
            
            fprintf('\nStarting to decode...\n');
            % Decode
            [uhat,xhat] = gf_fft_decode_qec(recv,f,p,m,[],lengths,N);
            uhat = uhat';
            fprintf('\nNumber of decoding errors = %d\n\n',nnz(u-uhat));
            biterr(i) = mean(uhat~=u);
        end
        % Display average bit and block error rate
        % mean(biterr)
        % mean(biterr>0)
        blockerr = [blockerr mean(biterr>0)];
    end
    fprintf('Code Rate: %f\n',rate);
    plot(0.1:0.1:1,blockerr);
    xlabel('Erasure Rate','FontSize',16);
    ylabel('Block Erasure Rate','FontSize',16);
    t = strcat('GF-FFT Code Performance; Blockerr = ',num2str(d));
    t = strcat(t,'; e = ',num2str(0.5),'; N = ',num2str(N),'; q = ',num2str(q),'; M = ',num2str(M),'; Rate = ',num2str(rate));
    title(t,'FontSize',16);
    t2 = clock;
    fprintf('Total Run Time: %f mins.\n',etime(t2,t1)/60);
else
    fprintf('\nN = %d does not divide p^m-1 = %d. Check.\n',N,p^m-1);
end