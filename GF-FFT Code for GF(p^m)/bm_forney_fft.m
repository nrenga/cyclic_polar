function [inp_est,r_est] = bm_forney_fft(r,inp,p,m)
% Run Berlekamp-Massey's algorithm to find lambda(x)
% Run Forney's algorithm to correct 'r' based on available data in
% 'inp'; 'inp' contains a few known values that should be obtained
%  when 'r' is evaluated at corresponding initial powers of
%  omega^-1, starting from 0 i.e. @ (omega^-1)^0,(omega^-1)^1 etc.

prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);
len_r = length(r);
q = p^m;

omega = 1;
for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
    % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
    a(1) = i;
    for j = 2:len_r
        a(j) = gfmul(a(j-1),i,field);
    end
    if (any(a==0) && find(a==0,1,'first') == len_r)
        omega = i;  % primitive len_r'th root of unity in GF(p^m); 'p' prime
        %omega = gfdiv(0,omega,field);   % Uncomment for RS codes
        %fprintf('%s%d\n','The primitive Nth root of unity was chosen to be alpha^',omega);
        break;
    end
end

erased_indices = find(r==1/q);  % Erased indices will have value = 1/q
r_est = r';
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
%fprintf('\nErasure Locator:');
%disp(tau);

no_of_syndromes = length(find(inp~=1/q));  % Number of known evaluations

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
    %fprintf('\nSyndromes:');
    %disp(S);
    
    if (~all(S == -Inf))    % Proceed if received word isn't a codeword
        % Find error locator polynomial
        phi = berlmass_fft(S,tau,p,m);
        %fprintf('\nError-Erasure Locator:');
        %disp(phi);
        phi_roots = find_roots_chien(phi,p,m);
        indices = gfdiv(zeros(1,length(phi_roots)),phi_roots,field)/omega + 1;
        
        % Error/erasure evaluator polynomial
        xt = [(-Inf)*ones(1,no_of_syndromes) 0];   % x^(2t)
        [quot, err_eval] = gfdeconv(gfconv(S,phi,field),xt,field);
        %fprintf('\nError-Erasure Evaluator:');
        %disp(err_eval);
        
        % Find Phi'(x)
        % phi = gfconv(lambda,tau,field);
        phi_d = gf_diff_poly(phi,field);
        
        omega_inv = gfdiv(0,omega,field);
        
        %fprintf('\nErasure Values:');
        % Compute error/erased values and correct r(x)
        for k = 1:length(indices)
            %fprintf('Index %d\n',indices(k));
            X_k_inv = gfpow(omega_inv,indices(k)-1,field);  % X_k^{-1} = (omega^{-i_k})
            err_eval_k = gf_eval_poly(err_eval,X_k_inv,field);
            phi_d_k = gf_eval_poly(phi_d,X_k_inv,field);
            
            err_val_k = gfsub(-Inf,gfdiv(err_eval_k,phi_d_k,field),field);
            if (b == 0)
                X_k = gfdiv(0,X_k_inv,field);
                err_val_k = gfmul(err_val_k,X_k,field);
            end
            r_est(indices(k)) = gfsub(r_est(indices(k)),err_val_k,field);
        end
        %fprintf('Decoded Codeword:');
        %disp(r_est');
    end
    % Compute the input based on corrected r(x)
    inp_est = gf_fft(r_est,p,m,[],length(r_est));
else
    inp_est = gf_fft(r_est,p,m,[],length(r_est));   % Will return garbage estimate
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