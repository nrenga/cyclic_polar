clc
close all
clear all
t1 = clock;
% Setup code parameters
n = 8;
N = 2^n;
% Find best frozen bits and compute rate
e = 0.5;
d = 0.1;   % Design Block Error Rate
f = polar_design_bec(n,e,d);
k = nnz(f==1/2);
rate = k/N;
individual_bit_errors = zeros(1,N);
blockerr = [];
for e = 0.45:0.1:0.65
    % Run a few sims to compare with union bound
    M = 1000;
    biterr = zeros(1,M);
    for i = 1:M
        clc
        fprintf('%s%f\n','Erasure Rate: ',e);
        fprintf('%s%d%s%d\n','Simulation ',i,' out of ',M);
        disp(blockerr);
        % Set frozen bits, add random data, and encode
        u = f;
        u(f==1/2) = rand(1,k)<0.5;
        x = polar_transform(u);
        % Transmit
        y = x;
        y(rand(1,N)<e)=1/2;  % In P1 domain, "0"=>P(1)=0, "1"=>P(1)=1, "1/2"=>P(1)=P(0)=1/2
        % Decode
        [uhat,xhat] = polar_decode(y,f);
        individual_bit_errors = individual_bit_errors + (uhat~=u);
        biterr(i) = mean(uhat~=u);
    end
    % Display average bit and block error rate
    % mean(biterr)
    blockerr = [blockerr mean(biterr>0)];
    individual_bit_errors = individual_bit_errors/M;
end
fprintf('Code Rate: %f\n',rate);
plot(0.45:0.1:0.65,blockerr);
xlabel('Erasure Rate','FontSize',16);
ylabel('Block Erasure Rate','FontSize',16);
t = strcat('Binary Polar Code Performance; Blockerr = ',num2str(d));
t = strcat(t,'; e = ',num2str(0.5),'; N = ',num2str(N),'; Rate = ',num2str(rate));
title(t,'FontSize',16);
t2 = clock;
fprintf('%s%f%s\n','Total Run Time: ',etime(t2,t1)/60,' mins.');