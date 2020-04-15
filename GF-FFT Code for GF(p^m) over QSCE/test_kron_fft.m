function [F,T,W,D,P] = test_kron_fft(f)
% function [F,T,W,D,P] = test_kron_fft(f)
%
% Test Cooley-Tukey Kronoecker Product Decomposition
%

% Check shuffle orientation (make sure code matches math in paper)
if 0
a = rand(2,2);
b = rand(3,3);
ab = kron(a,b);
ba = kron(b,a);
P = pmat(2,3);
P'*ab*P-ba
P*[1:6]'
pause
end

% Generate kronecker fft
n = length(f);
N = prod(f);
P = cell(1,n);
D = cell(1,n);
W = cell(1,n);
T = cell(1,n);
p = 1;
F = eye(N);
for i=1:n
  p = p*f(i);
  q = N/p;
  m = N/f(i);
  P{i} = pmat(q,f(i));
  D{i} = dmat(f(i),q);
  W{i} = dftmtx(f(i));
  T{i} = kron(P{i}*D{i},eye(p/f(i))) * kron(W{i},eye(m));
  F = T{i}*F;
end

% Generate P matrix
function p = pmat(a,b)

index = 0:(a*b-1);
index2 = reshape(reshape(index,[a b])',1,a*b);
p = sparse(index+1,index2+1,ones(1,a*b));
return

% Generate D matrix
function d = dmat(a,b)
  
index = 0:(a*b-1);
s = floor(index/b);
r = index-s*b;
d = diag(exp(-2*pi*1i/a/b).^(s.*r));
return
 
