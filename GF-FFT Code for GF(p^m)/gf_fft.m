function [ y ] = gf_fft(u,p,m,omega,N)
% FFT on a Galois Field GF(q=p^m)
% field contains elements of GF(q=p^m) in the exponential format
% N is the length of input vector which should be a power of 2
% Assumes N|q-1 and that q is prime
% The parameter omega is optional. The code anyway finds and uses the
% smallest primitive Nth root of unity in GF(q=p^m)
% All values are in exponential format with respect to the primitive
% element of GF(p^m)

len = length(u);  % A power of 2 that divides q-1
omega_valid = 1;
lengths = factor(len);
prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);

if (mod(p^m-1,N) == 0)
    if (N == len)
        omega_valid = 0;
        for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
            % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
            a(1) = i;
            for j = 2:N
                a(j) = gfmul(a(j-1),i,field);
            end
            if (any(a==0) && find(a==0,1,'first') == N)
                omega = i;
                %fprintf('%s%d\n','The primitive Nth root of unity was chosen to be alpha^',omega);
                omega_valid = 1;
                break;
            end
        end
    end
    if (omega_valid == 1)
        if (isprime(len))
            %disp(u);
            y = naive_gf_fft(u,p,m);
        else
            no_of_blocks = prod(lengths(2:end));
            u_inp = (reshape(u,lengths(1),no_of_blocks))';
            z = (-Inf)*ones(no_of_blocks,lengths(1));
            for i = 1:size(u_inp,2)
                ordered_inp = (u_inp(:,i))';
                z(:,i) = gf_fft(ordered_inp,p,m,omega,N);
            end
            % Multiply each element by omega^(i'k'')
            gamma = omega*(N/(size(z,1)*size(z,2)));
            for i = 0:(size(z,2)-1)
                for k = 0:(size(z,1)-1)
                    z(k+1,i+1) = gfmul(z(k+1,i+1),mod(gamma*i*k,p^m-1),field);
                end    
            end
            y = (-Inf)*ones(no_of_blocks,lengths(1));
            for i = 1:size(z,1)
                y(i,:) = (gf_fft(z(i,:),p,m,omega,N))';
            end      
            y = (reshape(y,1,[]))';
        end
    else
        fprintf('%s\n','Unable to find a primitive Nth root of unity in this field. GF-FFT Aborted');
    end
else
    fprintf('%s\n','The input vector length should divide (q-1)! GF-FFT Aborted');
end