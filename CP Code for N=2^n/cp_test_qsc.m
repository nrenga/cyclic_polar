clc
close all
clear all
t1 = clock;
% Setup code parameters
n = 4;
N = 2^n;
p = 0.5;
d = 0.1;   % Design Block Error Rate
q = 257;
f = cp_design_qec(n,p,d,q);  % Design is same with and without twiddle
k = nnz(f(1,:)~=1);
rate = k/N;

% Find best frozen bits and compute rate
% f = input a priori probs matrix in input order;
%     each column signifies a particular bit; each row signifies a GF value
M = 1000;
%blockerr = [0 0 0 0.12 0.29 0.38 0.55 0.73 0.89 1 1 1 1 1]; N=256
%blockerr = [0 0.01 0.02 0.06 0.11 0.15 0.2 0.24 0.28 0.58 0.91 0.94 1 1]; N=16
blockerr = [];
for p = [0.35]
    % Run a few sims
    biterr = zeros(1,M);
    for i = 1:M
        clc
        fprintf('Error Rate: %f\n',p);
        fprintf('Simulation %d out of %d\n',i,M);
        disp(blockerr);
        % Set frozen bits, add random data, and encode
        u = f(1,:);
        u(u==1/q) = randi(q,1,k)-1;
        u(u==1) = 0;
        x = gf_fft(u,q,[],N);
        % Transmit through q-ary symmetric channel QSC(p)
        
        recv = x;
        %recv(rand(1,N)<e) = 1/q;  % Random erasures
        error_indices = find(rand(1,N)<p);
        for j = 1:length(error_indices)
            set = 0:(q-1);
            set(set==x(error_indices(j))) = [];
            recv(error_indices(j)) = set(randi(q-1));
        end
        
        % y = bit APP matrix from channel in output order;
        %     each column signifies a particular bit; each row signifies a GF value
        y = zeros(q,N);
        for j = 1:N
            %if (recv(j) == 1/q)
                y(:,j) = (p/(q-1))*ones(q,1);
                y(recv(j)+1,j) = 1-p;    % Received value is correct with prob. (1-p)
            %else
            %    y(:,j) = zeros(q,1);
            %    y(recv(j)+1,j) = 1;  % Unerased bit; Set the probability for that value as 1
            %end
        end
        % Decode
        [uhat,xhatprob] = cp_decode(y,f,q,[],[],N);
        %         xhat = zeros(1,N);
        %         for j = 1:N
        %             if (recv(j) == 1/q)
        %                 ind = find(xhatprob(:,j) == max(xhatprob(:,j)),1,'first');
        %                 xhat(1,j) = ind - 1;
        %             else
        %                 xhat(1,j) = recv(j);
        %             end
        %         end
        biterr(i) = mean(uhat~=u);
    end
    % Display average bit and block error rate
    % mean(biterr)
    % mean(biterr>0)
    blockerr = [blockerr mean(biterr>0)];
end
fprintf('Code Rate: %f\n',rate);
p = [0.35];
%p = [0.1 0.2 0.3 0.4 0.42 0.44 0.46 0.48 0.5 0.6:0.1:1];
%p = 0.1:0.1:1;
plot(p,blockerr);
xlabel('Error Rate','FontSize',16);
ylabel('Block Error Rate','FontSize',16);
t = strcat('CP Code Performance; Blockerr = ',num2str(d));
t = strcat(t,'; e = ',num2str(0.5),'; N = ',num2str(N),'; q = ',num2str(q),'; M = ',num2str(M),'; Rate = ',num2str(rate));
title(t,'FontSize',16);
t2 = clock;
fprintf('Total Run Time: %f mins.\n',etime(t2,t1)/60);