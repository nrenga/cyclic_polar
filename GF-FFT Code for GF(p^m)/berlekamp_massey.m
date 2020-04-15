% Berlekamp-Massey Algorithm to decode RS and BCH Codes (Non-binary)

clc
close all
clear all
p = 2;
fprintf('%s\n\n','This code runs the B-M Algorithm only for GF(2^m)');
% n = input('Give the length n = 2^m - 1 of the code: ');
% n = 15;
n = 255;
m = log2(n+1);
fprintf('%s%d%s\n','The field we are working on is GF(2^',m,')');
prim_poly = gfprimdf(m);
fprintf('%s%s\n','The default primitive polynomial to generate this field is ',convert2poly(prim_poly));
% S = input('Give the exponents of the syndromes in a row matrix: ');
% S = [13 11 0 7 0 0];
S = [77 102 170 3 159 171 122 189 16 32 151 75 241 192 165 154];
N = length(S);
t = N/2;
field = gftuple([-1:p^m-2]',prim_poly,p);

%Algorithm Begins
L = 0;
cx = 0; % '0' implies the element '1'; Connection Polynomial
px = 0; % '0' implies the element '1'; Connection Polynomial before last length change
l = 1; % l is (k-m), the amount of shift in update
dm = 0; %Previous discrepancy

for k = 1:N
    d =  S(1,k);
    if (length(cx) < (L+1))
        cx = [cx (-1)*ones(1,L+1-length(cx))];
    end
    for i = 1:L
        d = gfadd(d,gfmul(cx(1,i+1),S(1,k-i),field),field);
    end
    if (d == 0)
        l = l + 1;
    else
        if ((2*L) >= k)
            if (dm == -Inf)
                temp = -1;
            else
                temp = gfmul(d,gfdiv(0,dm,field),field);  % gfdiv(0,dm,field) = dm^(-1)
            end    
            xp = gfconv([(-1)*ones(1,l) temp],px,field);   % d*(dm^(-1))*(x^l)*p(x)
            if (length(cx) < length(xp))  % The length of matrix representing (x^l) is (l+1)
                cx = [cx (-1)*ones(1,length(xp)-length(cx))];
            end
            cx = gfsub(cx,xp,field);
            l = l + 1;
        else
            tx = cx;
            if (dm == -Inf)
                temp = -1;
            else
                temp = gfmul(d,gfdiv(0,dm,field),field);  % gfdiv(0,dm,field) = dm^(-1)
            end    
            xp = gfconv([(-1)*ones(1,l) temp],px,field);   % d*(dm^(-1))*(x^l)*p(x)
            if (length(cx) < length(xp))  % The length of matrix representing (x^l) is (l+1)
                cx = [cx (-1)*ones(1,length(xp)-length(cx))];
            end
            cx = gfsub(cx,xp,field);
            L = k - L;
            px = tx;
            dm = d;
            l = 1;
        end
    end
    fprintf('%s%d%s%d%s%d%s%d%s%d\n','k = ',k,'; dk = ',d,'; L = ',L,'; l = ',l,'; dm = ',dm);
%     fprintf('%s','p(x) = ');
%     disp(px);
    fprintf('\n%s\n\n','c(x) = ');
    disp(cx);
end