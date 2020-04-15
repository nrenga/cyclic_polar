function [lambda] = berlmass_fft(S,tau,p,m)
% Berlekamp-Massey Algorithm to decode FFT Codes (Non-binary)

prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);

N = length(S);

%Algorithm Begins
no_of_erasures = length(tau)-1;
L = length(tau) - 1;  % Current length of LFSR
cx = tau; % '0' implies the element '1'; Connection Polynomial
px = tau; % '0' implies the element '1'; Connection Polynomial before last length change
l = 1;  % l is (k-m), the amount of shift in update
dm = 0; % Previous discrepancy

for k = length(tau):N
    d =  S(1,k);
    if (length(cx) < (L+1))
        cx = [cx (-Inf)*ones(1,L+1-length(cx))];
    end
    for i = 1:L
        d = gfadd(d,gfmul(cx(1,i+1),S(1,k-i),field),field);
    end
    if (d == -Inf)
        l = l + 1;
    else
        if ((2*L) >= k + no_of_erasures)
            if (dm == -Inf)
                temp = -Inf;
            else
                temp = gfmul(d,gfdiv(0,dm,field),field);  % gfdiv(0,dm,field) = dm^(-1)
            end
            xp = gfconv([(-Inf)*ones(1,l) temp],px,field);   % d*(dm^(-1))*(x^l)*p(x)
            if (length(cx) < length(xp))  % The length of matrix representing (x^l) is (l+1)
                cx = [cx (-Inf)*ones(1,length(xp)-length(cx))];
            end
            cx = gfsub(cx,xp,field);
            l = l + 1;
        else
            tx = cx;
            if (dm == -Inf)
                temp = -Inf;
            else
                temp = gfmul(d,gfdiv(0,dm,field),field);  % gfdiv(0,dm,field) = dm^(-1)
            end
            xp = gfconv([(-Inf)*ones(1,l) temp],px,field);   % d*(dm^(-1))*(x^l)*p(x)
            if (length(cx) < length(xp))  % The length of matrix representing (x^l) is (l+1)
                cx = [cx (-Inf)*ones(1,length(xp)-length(cx))];
            end
            cx = gfsub(cx,xp,field);
            L = k - (L - no_of_erasures);
            px = tx;
            dm = d;
            l = 1;
        end
    end
    %fprintf('%s%d%s%d%s%d%s%d%s%d\n','k = ',k,'; dk = ',d,'; L = ',L,'; l = ',l,'; dm = ',dm);
    %     fprintf('%s','p(x) = ');
    %     disp(px);
    %fprintf('\n%s\n\n','c(x) = ');
    %disp(cx);
end
lambda = cx;
end