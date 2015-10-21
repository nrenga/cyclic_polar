clc
close all
clear all
t1 = clock;
% Setup code parameters
n = 8;
N = 2^n;
rate = 0.5;   % Design Rate of the Code
% Find best frozen bits and compute rate
k = rate*N;   % For good bits, "1" and "0" are equi-probable. So P(1) = 1/2
f = polar_rate_design_bec(n,k);
% e = 0.5;
% d = 0.1;
blockerr = [];
for e = 0.1:0.1:1
    % Run a few sims to compare with union bound
    M = 10;
    biterr = zeros(1,M);
    for i=1:M
        clc
        fprintf('%s%f\n','Erasure Rate: ',e);
        fprintf('%s%d%s%d\n','Simulation ',i,' out of ',M);
        % Set frozen bits, add random data, and encode
        u = f;
        u(f==1/2) = rand(1,k)<0.5;
        x = polar_transform(u);
        % Transmit
        y = x;
        y(rand(1,N)<e)=1/2;  % In P1 domain, "0"=>P(1)=0, "1"=>P(1)=1, "1/2"=>P(1)=P(0)=1/2
        % Decode
        [uhat,xhat] = polar_decode(y,f);
        biterr(i) = mean(uhat~=u);
    end
    % Display average bit and block error rate
    % mean(biterr)
    blockerr = [blockerr mean(biterr>0)];
end
plot(0.1:0.1:1,blockerr);
xlabel('Erasure Rate','FontSize',16);
ylabel('Block Erasure Rate','FontSize',16);
t = strcat('Binary Polar Code Performance; Design Rate = ',num2str(rate));
t = strcat(t,'; e = ',num2str(0.5),'; N = ',num2str(N));
title(t,'FontSize',16);
t2 = clock;
fprintf('%s%f%s\n','Total Run Time: ',etime(t2,t1)/60,' mins.');