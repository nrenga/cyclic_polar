% Code to test the CP code of blocklength N=2^n over QEC(e)
% This implementation has the message vector u in the time domain

clc
close all
clear all
t1 = clock;
% Setup code parameters
n = 4;
N = 2^n;
e = 0.5;
d = 0.1;   % Design Block Error Rate
q = 257;
f = cp_design_qec(n,e,d,q);  % Design is same with and without twiddle
k = nnz(f(1,:)~=1);
rate = k/N;

% Find best frozen bits and compute rate
% f = input a priori probs matrix in input order;
%     each column signifies a particular bit; each row signifies a GF value
M = 1000;
%blockerr = [0 0 0 0 0.09 0.74 1];
blockerr = [];
for e = 0.15
    % Run a few sims
    biterr = zeros(1,M);
    for i = 1:M
        clc
        fprintf('Erasure Rate: %f\n',e);
        fprintf('Simulation %d out of %d\n',i,M);
        disp(blockerr);
        % Set frozen bits, add random data, and encode
        u = f(1,:);
        u(u==1/q) = randi(q,1,k)-1;
        u(u==1) = 0;
        x = gf_fft(u,q,[],N);
        % Transmit
        
        recv = x;
        recv(rand(1,N)<e) = 1/q;  % Random erasures
        
        % y = bit APP matrix from channel in output order;
        %     each column signifies a particular bit; each row signifies a GF value
        y = zeros(q,N);
        for j = 1:N
            if (recv(j) == 1/q)
                y(:,j) = (1/q)*ones(q,1);
            else
                y(:,j) = zeros(q,1);
                y(recv(j)+1,j) = 1;  % Unerased bit; Set the probability for that value as 1
            end
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
plot(0.65,blockerr);
xlabel('Erasure Rate','FontSize',16);
ylabel('Block Erasure Rate','FontSize',16);
t = strcat('CP Code Performance; Blockerr = ',num2str(d));
t = strcat(t,'; e = ',num2str(0.5),'; N = ',num2str(N),'; q = ',num2str(q),'; M = ',num2str(M),'; Rate = ',num2str(rate));
title(t,'FontSize',16);
t2 = clock;
fprintf('Total Run Time: %f mins.\n',etime(t2,t1)/60);