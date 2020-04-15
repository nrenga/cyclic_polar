function [inp_est,r_est] = forney_fft(r,inp,p,m)
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
        %fprintf('%s%d\n','The primitive Nth root of unity was chosen to be alpha^',omega);
        break;
    end
end

omega_inv = gfdiv(0,omega,field);
erased_indices = find(r==1/q);  % Erased indices will have value = 1/q
r_est = r';
% b is the smallest exponent of omega_inv that is a root of g(x)
b = find(inp~=1/q,1,'first')-1;

if (~isempty(erased_indices))
    % Compute erasure locator polynomial
    lambda = 0;
    for l = 1:length(erased_indices)
        X_l = gfpow(omega_inv,erased_indices(l)-1,field);  % Xl = (omega^-1)^i_l
        fact = [0 gfsub(-1,X_l,field)];     % fact = (1-X_l*x)
        lambda = gfconv(lambda,fact,field); % lambda(x)=PI((1-X_l*x))
    end
    r_est(erased_indices) = -1;         % 0 in exponential format
    %fprintf('Received Word:');
    %disp(r);
    
    no_of_syndromes = length(find(inp~=1/q));  % Number of known evaluations
    
    if (no_of_syndromes > 0)
        % Find len_r^-1 in GF(p^m)
        n = gf_find_n(len_r,field);
        len_r_inv = gfdiv(0,n,field);
        
        S = (-1)*ones(1,max((b-1),0)+no_of_syndromes);   % Syndrome Polynomial
        for ns = b:(b+no_of_syndromes-1)
            point = gfpow(omega_inv,ns,field);
            r_eval = gfmul(len_r_inv,gf_eval_poly(r_est,point,field),field);
            if (b > 0)
                S(ns) = gfsub(r_eval,inp(ns+1),field);   % Syndrome
            else
                S(ns+1) = gfsub(r_eval,inp(ns+1),field);   % Syndrome
            end
        end
        %fprintf('\nSyndromes:');
        %disp(S);
        
        % Error evaluator polynomial
        xt = [(-1)*ones(1,no_of_syndromes) 0];   % x^(2t)
        [quot, err_eval] = gfdeconv(gfconv(S,lambda,field),xt,field);
        
        % Find lamba'(x)
        lambda_d = gf_diff_poly(lambda,field);
        
        %fprintf('\nErasure Values:');
        % Compute error values and correct r(x)
        for k = 1:length(erased_indices)
            X_k_inv = gfpow(omega,erased_indices(k)-1,field);  % Xl^-1 = (omega^i_l)
            err_eval_k = gf_eval_poly(err_eval,X_k_inv,field);
            lambda_d_k = gf_eval_poly(lambda_d,X_k_inv,field);
            err_val_k = gfsub(-1,gfdiv(err_eval_k,lambda_d_k,field),field);
            if (b == 0)
                X_k = gfdiv(0,X_k_inv,field);
                err_val_k = gfmul(err_val_k,X_k,field);
            end
            actual_err = gfmul(n,err_val_k,field);
            %disp(actual_err);
            r_est(erased_indices(k)) = gfsub(r_est(erased_indices(k)),actual_err,field);
        end
        %fprintf('Decoded Codeword:');
        %disp(r_est);

        % Compute the input based on corrected r(x)
        inp_est = naive_gf_ifft(r_est,p,m);
    else
        inp_est = naive_gf_ifft(r_est,p,m);   % Will return garbage estimate
    end
else
    inp_est = naive_gf_ifft(r_est,p,m);   % Exact inputs since no erasures
end
r_est(r_est == -Inf) = -1;
inp_est(inp_est == -Inf) = -1;
return
end

    function [val] = gf_eval_poly(px,point,field)
        % Evaluate the polynomial px at the value 'point' in GF(p^m)
        % 'point' is in exponential format
        
        val = -1;
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
        res = -1;       % -1 is 0 in exponential format
        for times = 1:n
            res = gfadd(res,0,field);  % Compute n in GF(p^m)
        end
        return
    end