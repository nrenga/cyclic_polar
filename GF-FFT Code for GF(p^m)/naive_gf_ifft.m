function [ y ] = naive_gf_ifft( u,p,m )
% O(N^2) naive implementation of IFFT just for verfiying the O(NlogN) IFFT
% IFFT on a Galois Field GF(p^m)
% N is the length of input vector which should be a power of 2
% Assumes N|p^m-1 and that p is prime
% All values are in exponential format with respect to the primitive
% element of GF(p^m)

N = length(u);
if (mod(p^m-1,N) == 0)
    prim_poly =gfprimfd(m,'min',p);
    field = gftuple([-1:(p^m-2)]',prim_poly,p);
    
    alpha = 1;   % Exponential format
    alpha_valid = 0;
    for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
        % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
        a(1) = i;
        for j = 2:N
            a(j) = gfmul(a(j-1),i,field);
        end
        if (any(a==0) && find(a==0,1,'first') == N)
            alpha = i;
            alphainv = 0;
            for j = 1:(p^m-2)
                alphainv = gfmul(alphainv,alpha,field);
            end
            alpha = alphainv;
            alpha_valid = 1;
            %fprintf('%s%d\n','The primitive Nth root of unity was chosen to be alpha^-',alpha);
            break;
        end
    end
    if (alpha_valid == 1)
        W = (-Inf)*ones(N,N);
        for i = 1:N
            % W(i,:) = mod((alpha^(i-1)).^[0:(N-1)],q);
            W(i,1) = 0;    % 0 is 1 in exponential format
            for j = 2:N
                % alpha1 = alpha^(i-1);
                alpha1 = 0;
                for k = 1:(i-1)
                    alpha1 = gfmul(alpha1,alpha,field);
                end
                W(i,j) = gfmul(W(i,j-1),alpha1,field);
            end
        end
        % Find N^-1 in GF(p^m)
        n = -Inf;       % -1 is 0 in exponential format
        for i = 1:N
            n = gfadd(n,0,field);  % Compute N in GF(p^m)
        end
        Ninv = gfdiv(0,n,field);
        
        % y = (mod(W*u',q))';
        y = (-Inf)*ones(N,1);   % -1 is 0 in exponential format
        for i = 1:N
            for j = 1:N
                y(i,1) = gfadd(y(i,1),gfmul(W(i,j),u(j),field),field);
            end
            y(i,1) = gfmul(Ninv,y(i,1),field);
        end
    else
        fprintf('%s\n','Unable to find a primitive Nth root of unity in this field');
        y = [];
    end
else
    fprintf('%s\n','N should divide p^m-1');
    y = [];
end
end