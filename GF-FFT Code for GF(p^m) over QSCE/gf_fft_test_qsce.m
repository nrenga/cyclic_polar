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
verbose = 0;

if (mod(q-1,N) == 0)
    lengths = factor(N);
    f = gf_fft_design_qec(lengths,e,d,q);  % Design is for the QEC
    k = nnz(f == 1/q);
    rate = k/N;
    
    % Find best frozen bits and compute rate
    % f = input a priori probs matrix in input order;
    %     each column signifies a particular bit; each row signifies a GF value
    M = 100;
    blockerr = [];
    %blockerr = [0 0.02 0.03 0.2 0.44 0.61 0.89 0.99 1 1 1 1]; % N = 63 over GF(64) on QSC e = 0.1:0.05:0.7
    %blockerr = [0 0 0 0 0.05 0.35 0.86 1 1 1]; % N = 63 over GF(64) on QEC e = 0.1:0.1:1
    %blockerr = [0 0.06 0.07 0.27 0.5 0.63 0.9 1 1]; % N = 255 over GF(256) on QSC e = [0.1 0.2:0.02:0.3 0.4 0.5]
    e_test = 0.7:0.1:1;
    for e = e_test
        % Run a few sims
        biterr = zeros(1,M);
        for i = 1:M
            clc
            fprintf('\nError Rate: %f\n',e);
            fprintf('\nSimulation %d out of %d\n',i,M);
            fprintf('\nBlock Error Rate:\n');
            disp(blockerr);
            
            % Set frozen bits, add random data, and encode
            u = (-Inf)*ones(1,N);  % In exponential format;
            u(f == 1/q) = randi(q,1,k)-2;
            u(u == -1) = -Inf;
            
            fprintf('\nEncoding now...\n');
            v = gf_ifft(u,p,m,[],N);
            
            r = v;
            
            % QSCE Channel
            perror = e;
            perasure = 0;
            pcorrect = 1 - perror - perasure;
            symbols = [1/q -Inf 0:(q-2)];
            for j = 1:N
                pvec = [perasure (perror/(q-1))*ones(1,q)];
                pvec(symbols == v(j)) = pcorrect;
                plen = length(pvec);     % plen = q+1
                pvec = reshape(pvec,plen,1);   % pvec = pvec'
                index = sum(repmat(rand(1),plen,1) > cumsum(pvec),1) + 1;
                r(j) = symbols(index);
            end

            %r(rand(1,N)<e) = 1/q;  % Random erasures
            
            erased_indices = find(r == 1/q);
            fprintf('\nNumber of erasures in received word = %d\n',nnz(erased_indices));
            
            unerased_indices = 1:N;
            unerased_indices(erased_indices) = [];
            
            error_indices = unerased_indices(r(unerased_indices) ~= v(unerased_indices));
            fprintf('\nNumber of errors in received word = %d\n',nnz(error_indices));
            
            fprintf('\nStarting to decode...\n');
            % Now decode
            [u_est,r_est] = gf_fft_decode_qsce(r,f,p,m,[],lengths,N,verbose);
            u_est = u_est';
          
            fprintf('\nNumber of decoding errors = %d\n\n',nnz(r_est ~= v));
            biterr(i) = mean(r_est ~= v);
        end
        % Display average bit and block error rate
        % mean(biterr)
        % mean(biterr>0)
        blockerr = [blockerr mean(biterr>0)];
    end
    fprintf('Code Rate: %f\n',rate);
    plot(e_test,blockerr);
    xlabel('Error Rate','FontSize',16);
    ylabel('Block Error Rate','FontSize',16);
    t = strcat('Cyclic Polar Code Performance on QSC; Blockerr = ',num2str(d));
    t = strcat(t,'; e = ',num2str(0.5),'; N = ',num2str(N),'; q = ',num2str(q),'; M = ',num2str(M),'; Rate = ',num2str(rate));
    title(t,'FontSize',16);
    t2 = clock;
    fprintf('Total Run Time: %f mins.\n',etime(t2,t1)/60);
else
    fprintf('\nN = %d does not divide p^m-1 = %d. Check.\n',N,p^m-1);
end