function [inp_est,r_est] = bm_forney_fft(r,inp,p,m,verbose)
% Run Berlekamp-Massey's algorithm to find lambda(x)
% Run Forney's algorithm to correct 'r' based on available data in
% 'inp'; 'inp' contains a few known values that should be obtained
%  when 'r' is evaluated at corresponding initial powers of
%  omega, starting from 0 i.e. @ omega^0,omega^1 etc.

prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);
N = length(r);

omega = 1;
for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
    % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
    a(1) = i;
    for j = 2:N
        a(j) = gfmul(a(j-1),i,field);
    end
    if (any(a==0) && find(a==0,1,'first') == N)
        omega = i;  % primitive Nth root of unity in GF(p^m); 'p' prime
        %fprintf('\nThe primitive root was chosen to be alpha^%d\n',omega);
        break;
    end
end

if (verbose >= 1)
    %disp(omega);
    fprintf('\nr = ');
    disp(r');
    fprintf('\ninp = ');
    disp(inp');
end
q = p^m;
erased_indices = find(r==1/q);  % Erased indices will have value = 1/q
r_est = r;
% b is the smallest exponent of omega_inv that is a root of g(x)
b = find(inp~=1/q,1,'first')-1;

% Compute erasure locator polynomial
tau = 0;
for l = 1:length(erased_indices)
    X_l = gfpow(omega,erased_indices(l)-1,field);  % Xl = omega^i_l
    fact = [0 gfsub(-Inf,X_l,field)];     % fact = (1-X_l*x)
    tau = gfconv(tau,fact,field);       % Tau(x) = PI((1-X_l*x))
end
r_est(erased_indices) = -Inf;         % 0 in exponential format
if (verbose >= 1)
    fprintf('\nErasure Locator:');
    disp(tau);
end

no_of_syndromes = length(find(inp~=1/q));  % Number of known evaluations

no_of_errors = 0;
no_of_erasures = length(erased_indices);

decoder_failure = 1;

if (no_of_syndromes > 0)
    S = (-Inf)*ones(1,max((b-1),0)+no_of_syndromes);   % Syndrome Polynomial
    for ns = b:(b+no_of_syndromes-1)
        point = gfpow(omega,ns,field);
        r_eval = gf_eval_poly(r_est,point,field);
        if (b > 0)
            S(ns) = gfsub(r_eval,inp(ns+1),field);   % Syndrome
        else
            S(ns+1) = gfsub(r_eval,inp(ns+1),field);   % Syndrome
        end
    end
    if (verbose == 2)
        fprintf('\nSyndromes:');
        disp(S);
    end
    
    if (~all(S == -Inf))    % Proceed if received word isn't a codeword
        % Find error locator polynomial
        phi = berlmass_fft(S,tau,p,m);
        if (verbose >= 1)
            fprintf('\nError-Erasure Locator:');
            disp(phi);
        end
        phi_roots = find_roots_chien(phi,p,m);
        
        no_of_erasures = length(erased_indices);
        no_of_errors = length(phi_roots) - no_of_erasures;
        
        if (~isempty(phi_roots) && no_of_syndromes >= 2*no_of_errors + no_of_erasures)
            
            indices = gfdiv(zeros(1,length(phi_roots)),phi_roots,field);
            
            if (all(mod(indices,omega) == 0))
                
                indices = indices/omega + 1;
                
                % Error/erasure evaluator polynomial
                xt = [(-Inf)*ones(1,no_of_syndromes) 0];   % x^(2t)
                [~, err_eval] = gfdeconv(gfconv(S,phi,field),xt,field);
                if (verbose == 2)
                    fprintf('\nError-Erasure Evaluator:');
                    disp(err_eval);
                end
                
                % Find Phi'(x)
                % phi = gfconv(lambda,tau,field);
                phi_d = gf_diff_poly(phi,field);
                
                omega_inv = gfdiv(0,omega,field);
                
                decoder_failure = 0;
                %fprintf('\nErasure Values:');
                % Compute error/erased values and correct r(x)
                for k = 1:length(indices)
                    %fprintf('Index %d\n',indices(k));
                    X_k_inv = gfpow(omega_inv,indices(k)-1,field);  % X_k^{-1} = (omega^{-i_k})
                    err_eval_k = gf_eval_poly(err_eval,X_k_inv,field);
                    phi_d_k = gf_eval_poly(phi_d,X_k_inv,field);
                    
                    if (phi_d_k == -Inf)
                        decoder_failure = 1;
                    else
                        err_val_k = gfsub(-Inf,gfdiv(err_eval_k,phi_d_k,field),field);
                        if (b == 0)
                            X_k = gfdiv(0,X_k_inv,field);
                            err_val_k = gfmul(err_val_k,X_k,field);
                        end
                        r_est(indices(k)) = gfsub(r_est(indices(k)),err_val_k,field);
                    end
                end
                
            end
        end
    else
        decoder_failure = 0;  % Received word is a codeword
    end
else
    % No input information is available. But, received word might be a
    % codeword. So, just give out an estimate based on it. If it had
    % errors, we are anyway going to make a mistake; so this estimate
    % does not harm the decoder performance.
    decoder_failure = 0;   % Return some estimate; might not be correct
end
% Compute the input based on corrected r(x)
if (decoder_failure == 0)
    if (no_of_syndromes == 0 && no_of_erasures > 0)
        r_est(erased_indices) = 1/q;
        inp_est = [inp; (1/q)*ones(length(r_est)-length(inp),1)];
    else
        inp_est = gf_fft(r_est,p,m,[],length(r_est));
        inp_est(inp ~= 1/q) = inp(inp ~= 1/q);   % Just to be safe
    end
else
    r_est(erased_indices) = 1/q;
    inp_est = [inp; (1/q)*ones(length(r_est)-length(inp),1)];
end
if (verbose >= 1)
    fprintf('Decoded Codeword:');
    disp(r_est');
    fprintf('Decoded Input:');
    disp(inp_est');
    fprintf('\n--------------------------------');
end
return
end

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
            fprintf('\nif case...\n');
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
            fprintf('\nelse case...\n');
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
    fprintf('%s%d%s%d%s%d%s%d%s%d\n','k = ',k,'; dk = ',d,'; L = ',L,'; l = ',l,'; dm = ',dm);
    fprintf('%s','p(x) = ');
    disp(px);
    fprintf('\n%s\n\n','c(x) = ');
    disp(cx);
end
lambda = cx;
end

function [poly_roots] = find_roots_chien(poly,p,m)
% Find the roots of the polynomial 'poly' in GF(p^m)

prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);

poly_roots = [];
no_of_roots = 0;
for point = [-Inf 0:(p^m-2)]
    %gf_eval_poly(poly,point,field)
    if (gf_eval_poly(poly,point,field) == -Inf)
        no_of_roots = no_of_roots + 1;
        poly_roots = [poly_roots point];
    end
end
return
end

function [val] = gf_eval_poly(px,point,field)
% Evaluate the polynomial px at the value 'point' in GF(p^m)
% 'point' is in exponential format

val = -Inf;
for k = 0:(length(px)-1)
    multiplier = gfpow(point,k,field);
    val = gfadd(val,gfmul(px(k+1),multiplier,field),field);
end
return
end

function [px_d] = gf_diff_poly(px,field)
% Differentiate the polynomial px in GF(p^m)

px_d = px(2:end);
for k = 1:length(px_d)
    multiplier = gf_find_n(k,field);
    px_d(k) = gfmul(px_d(k),multiplier,field);
end
return
end

function [res] = gfpow(a,b,field)
% Computes a^b in GF(p^m) represented as field
res = 0;   % 0 is 1 in exponential format
for times = 1:b
    res = gfmul(res,a,field);
end
return
end

function [res] = gf_find_n(n,field)
% Find the value of 'n' in GF(p^m)
res = -Inf;       % -1 is 0 in exponential format
for times = 1:n
    res = gfadd(res,0,field);  % Compute n in GF(p^m)
end
return
end