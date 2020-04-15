function [ x ] = naive_gf_fft( u,q )
% O(N^2) naive implementation of FFT just for verfiying the O(NlogN) FFT
% FFT on a Galois Field GF(q)
% N is the length of input vector which should be a power of 2
% Assumes N|q-1 and that q is prime

N = length(u);
alpha = 2;   % Assume 2 is primitive root and then check
alpha_valid = 0;
for i = 2:(q-1)   % Find the smallest primitive Nth root of unity
    % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
    a(1) = i;
    for j = 2:N
        a(j) = mod(a(j-1)*i,q);
    end
    if (any(a==1) && find(a==1,1,'first') == N)
        alpha = i;
        alpha_valid = 1;
        fprintf('%s%d\n','The primitive root was chosen to be ',alpha);
        break;
    end
end
if (alpha_valid == 1)
    W = zeros(N,N);
    for i = 1:N
        % W(i,:) = mod((alpha^(i-1)).^[0:(N-1)],q);
        W(i,1) = 1;
        for j = 2:N
            % alpha1 = alpha^(i-1);
            alpha1 = 1;
            for k = 1:(i-1)
                alpha1 = mod(alpha1*alpha,q);
            end
            W(i,j) = mod(W(i,j-1)*alpha1,q);
        end
    end
    if (size(u,2) == 1)
        u = u';
    end
    % x = (mod(W*u',q))';
    for i = 1:N
        x(1,i) = 0;
        for j = 1:N
            x(1,i) = mod(x(1,i)+mod(W(i,j)*u(j),q),q);
        end
    end
else
    fprintf('%s\n','Unable to find a primitive Nth root of unity in this field');
end
end