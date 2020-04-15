function [ y ] = naive_gf_fft( u,p,m )
% O(N^2) naive implementation of FFT just for verfiying the O(NlogN) FFT
% FFT on a Galois Field GF(q=p^m)
% N is the length of input vector which should be a power of 2
% Assumes N|p^m-1 and that p is prime
% All values are in exponential format with respect to the primitive
% element of GF(p^m)

N = length(u);
prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);

alpha = 1;   % exponential format
alpha_valid = 0;
for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
    % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
    a(1) = i;
    for j = 2:N
        a(j) = gfmul(a(j-1),i,field);
    end
    if (any(a==0) && find(a==0,1,'first') == N)
        alpha = i;
        alpha_valid = 1;
        %fprintf('%s%d\n','The primitive Nth root of unity was chosen to be alpha^',alpha);
        break;
    end
end
if (alpha_valid == 1)
    W = zeros(N,N);
    for i = 1:N
        % W(i,:) = mod((alpha^(i-1)).^[0:(N-1)],q);
        W(i,1) = 0;     % exponential format; W(i,1) = alpha^0 = 1
        for j = 2:N
            % alpha1 = alpha^(i-1);
            alpha1 = 0;     % exponential format
            for k = 1:(i-1)
                alpha1 = gfmul(alpha1,alpha,field);
            end
            W(i,j) = gfmul(W(i,j-1),alpha1,field);
        end
    end
    if (size(u,2) == 1)
        u = u';
    end
    % y = (mod(W*u',q))'
    y = (-Inf)*ones(N,1);    % Y(x,x) = alpha^-1 = 0
    for i = 1:N
        for j = 1:N
            y(i,1) = gfadd(y(i,1),gfmul(W(i,j),u(j),field),field);
        end
    end
else
    fprintf('%s\n','Unable to find a primitive Nth root of unity in this field');
end
end