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

    function [res] = gfpow(a,b,field)
        % Computes a^b in GF(p^m) represented as field
        res = 0;   % 0 is 1 in exponential format
        for times = 1:b
            res = gfmul(res,a,field);
        end
        return
    end
    