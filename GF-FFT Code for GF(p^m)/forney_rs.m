function [rc] = forney_rs(r,t,lambda,roots,p,m)
% Run Forney's algorithm to correct 'r'

len_r = length(r);
prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);

omega = 1;
for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
    % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
    a(1) = i;
    for j = 2:len_r
        a(j) = gfmul(a(j-1),i,field);
    end
    if (any(a==0) && find(a==0,1,'first') == len_r)
        omega = i;  % primitive len_r'th root of unity in GF(p^m); 'p' prime
        fprintf('%s%d\n','The primitive root was chosen to be ',omega);
        break;
    end
end

no_of_syndromes = 2*t;
no_of_errors = length(lambda)-1;

S = (-1)*ones(1,no_of_syndromes);   % Syndrome Polynomial
for ns = 1:no_of_syndromes
    point = gfpow(omega,ns,field);
    S(ns) = gf_eval_poly(r,point,field);
end

% Error evaluator polynomial
xt = [(-1)*ones(1,no_of_syndromes) 0];   % x^(2t)
[quot, err_eval] = gfdeconv(gfconv(S,lambda,field),xt,field);

% Find lambda'(x)
lambda_d = gf_diff_poly(lambda,field);
%roots = gfroots(lambda,m,p);

% Compute error values and correct r(x)
for k = 1:no_of_errors
    X_l_inv = roots(k);  % Xl^-1 = (omega^-1)^i_l
    err_eval_k = gf_eval_poly(err_eval,X_l_inv,field);
    lambda_d_k = gf_eval_poly(lambda_d,X_l_inv,field);
    err_val_k = gfsub(-1,gfdiv(err_eval_k,lambda_d_k,field),field);
    err_loc_k = gfdiv(0,roots(k),field);
    r(err_loc_k+1) = gfsub(r(err_loc_k+1),err_val_k,field);
end
rc = r;
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
        % Find N^-1 in GF(p^m)
        res = -1;       % -1 is 0 in exponential format
        for times = 1:n
            res = gfadd(res,0,field);  % Compute n in GF(p^m)
        end
        return
    end