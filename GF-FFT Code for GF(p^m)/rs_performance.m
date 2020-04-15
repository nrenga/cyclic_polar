clc
close all
clear all

% (256,84) RS Code on QSC
k = 84;
N = 256;
t = floor((N-k)/2);

% QSC(e)
e = [0.1:0.05:0.35 0.4:0.02:0.5 0.6:0.1:1];
L = length(e);

blockerr = zeros(1,L);

for i = (t+1):N
    blockerr = blockerr + nchoosek(N,i) * ((e.^i).*((1-e).^(N-i)));
end
legend1 = '(256,84) extended RS Code';

% (256,84) Cyclic Polar Code over GF(257) with soft decoding on QSC
e2 = [0.1 0.2 0.3 0.35 0.4:0.02:0.5 0.6:0.1:1];
blockerr2 = [0 0 0.005 0.027 0.146 0.257 0.393 0.574 0.735 0.872 0.999 1 1 1 1];
legend2 = '(256,84) Cyclic Polar Code';

% (255,98) RS Code on QSC
k = 98;
N = 255;
t = floor((N-k)/2);

% QSC(e)
e3 = [0.1:0.05:0.35 0.4:0.02:0.5 0.6:0.1:1];
L = length(e3);

blockerr3 = zeros(1,L);

for i = (t+1):N
    blockerr3 = blockerr3 + nchoosek(N,i) * ((e3.^i).*((1-e3).^(N-i)));
end
legend3 = '(255,98) RS Code';

% (255,98) Cyclic Polar Code over GF(256) with hard decoding on QSC
e4 = [0.1 0.2:0.02:0.3 0.4 0.5 0.6:0.1:1];
blockerr4 = [0 0.06 0.07 0.27 0.5 0.63 0.9 1 1 ones(1,5)];
legend4 = '(255,98) Cyclic Polar Code';

% (16,4) RS Code on QSC
k = 4;
N = 16;
t = floor((N-k)/2);

% QSC(e)
e5 = [0.1:0.05:0.35 0.4:0.02:0.5 0.6:0.1:1];
L = length(e5);

blockerr5 = zeros(1,L);

for i = (t+1):N
    blockerr5 = blockerr5 + nchoosek(N,i) * ((e5.^i).*((1-e5).^(N-i)));
end
legend5 = '(16,4) RS Code';

% (255,98) Cyclic Polar Code over GF(256) with hard decoding on QSC
e6 = [0.1 0.2 0.3 0.35 0.4:0.02:0.5 0.6:0.1:1];
blockerr6 = [0 0.002 0.015 0.04 0.086 0.131 0.162 0.193 0.232 0.295 0.558 0.802 0.953 1 1];
legend6 = '(16,4) Cyclic Polar Code';

figure; 
semilogy(e,blockerr,'-or');
hold on
semilogy(e2,blockerr2,'-+r');
hold on
semilogy(e3,blockerr3,'-.ob');
hold on
semilogy(e4,blockerr4,'-.+b');
hold on
semilogy(e5,blockerr5,'--og');
hold on
semilogy(e6,blockerr6,'--xg');
ylim([0.001 1]);
legend(legend1, legend2, legend3, legend4, legend5, legend6);