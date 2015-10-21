function [ x ] = gf_ifft(u,q,omegainv,N)
% IFFT on a Galois Field GF(q)
% N is the length of input vector which should be a power of 2
% Assumes N|q-1 and that q is prime
% The parameter alpha is optional. The code anyway finds and uses the
% smallest primitive Nth root of unity in GF(q)

len = length(u);  % A power of 2 that divides q-1
omega_valid = 1;

if (mod(q-1,N) == 0)
    if (N == length(u))
        omega_valid = 0;
        for i = 2:(q-1)   % Find the smallest primitive Nth root of unity
            % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
            a(1) = i;
            for j = 2:N
                a(j) = mod(a(j-1)*i,q);
            end
            if (any(a==1) && find(a==1,1,'first') == N)
                omega = i;
                %fprintf('%s%d\n','The primitive root was chosen to be ',omega);
                omegainv = 1;
                for j = 1:(q-2)
                    omegainv = mod(omegainv*omega,q);  % omega^(q-1) = 1
                end
                %fprintf('%s%d\n','Thus, the inverse is ',omegainv);
                omega_valid = 1;
                break;
            end
        end
    end
    if (omega_valid == 1)
        if (len == 1)
            x = u;
        else
            %alpha = mod(omegainv^(N/length(u)),q);
            alpha = 1;
            for j = 1:(N/length(u))
                alpha = mod(alpha*omegainv,q);
            end
            % factors = mod(alpha.^[0:1:(len/2-1)],q);
            factors(1) = 1;
            for j = 2:(len/2)
                factors(j) = mod(factors(j-1)*alpha,q);
            end
            
            %Actual operation
            even = gf_ifft(u(1:2:end),q,omegainv,N);
            odd = gf_ifft(u(2:2:end),q,omegainv,N);
            x(1:len/2) = mod(even + mod(factors.*odd,q),q);
            % twiddle = mod(omegainv^(N/2),q);
            twiddle = 1;
            for j = 1:(N/2)
                twiddle = mod(twiddle*omegainv,q);
            end
            x((len/2+1):len) = mod(even + mod(twiddle*(mod(factors.*odd,q)),q),q);
            if (len == N)
                Ninv = gfdiv(1,N,q);
                x = mod(Ninv*x,q);
            end
        end
    else
        fprintf('%s\n','Unable to find a primitive Nth root of unity in this field');
    end
else
    fprintf('%s\n','The input vector length should divide (q-1)!');
end
end