function [ gx ] = gf_fft_find_gen( N,q,e,d )
% Find generator polynomial for the FFT polar code over the BEC
% N is the blocklength, code over GF(q) where q is prime and N|(q-1)
% e is the erasure rate of the BEC, d is the desired block error rate

f = gf_fft_design_bec(N,q,e,d);
indices = find(f==0)-1;

for i = 2:(q-1)   % Find the smallest primitive Nth root of unity
    % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
    a(1) = i;
    for k = 2:N
        a(k) = mod(a(k-1)*i,q);
    end
    if (any(a==1) && find(a==1,1,'first') == N)
        alphainv = find_gf_inv(i,q);
        %fprintf('%s%d\n','The inverse was found to be ',alphainv);
        break;
    end
end

gx = 1;
for i = 1:length(indices)
    factor = 1;
    for k = 1:indices(i)
        factor = mod(alphainv*factor,q);
    end
    gx = gfconv(gx,[mod((-1)*factor,q) 1],q);
end
return

% Find inverse of an element in GF(q)
    function ainv = find_gf_inv(a,q)
        ainv = 1;
        for j = 1:(q-1)
            if (mod(ainv*a,q) ~= 1)
                ainv = mod(ainv*a,q);
            else
                break;
            end
        end
        return
    end
end