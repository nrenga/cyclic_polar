function [u,x] = gf_fft_decode_qec(y,f,p,m,omega,lengths,N)
% y = received vector from channel in output order;
% f = input a priori probs vector in input decoding order;
%     each column signifies a particular bit; each row signifies a GF value
% u = input hard decisions in input order
% omega = primitive Nth root of unity in GF(p^m) where 'p' is prime

prim_poly =gfprimfd(m,'min',p);
field = gftuple([-1:(p^m-2)]',prim_poly,p);

if (N == prod(lengths))
    for i = 1:(p^m-2)   % Find the smallest primitive Nth root of unity
        % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
        a(1) = i;
        for j = 2:N
            a(j) = gfmul(a(j-1),i,field);
        end
        if (any(a==0) && find(a==0,1,'first') == N)
            omega = i;  % primitive Nth root of unity in GF(p^m); 'p' prime
            omega = gfdiv(0,omega,field);
            %fprintf('\nThe primitive root was chosen to be alpha^%d\n',omega);
            break;
        end
    end
end
% Recurse down to prime length
if (length(lengths) == 1)
    %fprintf('\nf = ');
    %disp(f');
    inp = f;    % Apriori knowledge of inputs
    if (all(inp == -Inf))   % All inputs are frozen
        u = inp;
        x = gf_ifft(u,p,m,[],length(u));
    else
        [u, x] = forney_fft(y,inp,p,m,mod(omega*(N/length(y)),p^m-1),field);
    end
    %fprintf('\ny = ');
    %disp(y);
    %fprintf('\nu = ');
    %disp(u');
    %fprintf('\nx = ');
    %disp(x');
    %fprintf('\n--------------------------------');
else
    len = lengths(1);
    no_of_blocks = prod(lengths(2:end));
    
    % Each column corresponds to the erasure probs. of a block
    y = (reshape(y,no_of_blocks,len))';
    x = (-Inf)*ones(len,no_of_blocks);       % input side of block
    
    for k = 1:no_of_blocks
        if (isempty(find(y(:,k) == 1/p^m,1)))   % Nothing erased
            x(:,k) = gf_fft(y(:,k),p,m,[],length(y(:,k)));
        else
            x(:,k) = 1/p^m;   % Pass back an erasure
        end
    end
    
    % Reorder input apriori information to input side FFT order
    old_f = f;
    f = [];
    for offset = 1:len
        f = [f old_f(offset:len:end)];
    end
    
    u = (-Inf)*ones(no_of_blocks,len);
    % Pass back to previous stage, decode and correct current stage
    for j = 1:len
        %disp(fft_order(:,j));
        % Multiply with gamma^(-i'*k'')
        gamma = mod(omega*(N/(size(x,1)*size(x,2))),p^m-1);  %omega*(N/(len*no_of_blocks))
        for k = 0:(size(x,2)-1)
            if (x(j,k+1) ~= 1/p^m)
                i = j-1;
                multiplier = gfdiv(0,mod(gamma*i*k,p^m-1),field);
                x(j,k+1) = gfmul(x(j,k+1),multiplier,field);
            end
        end
        
        [ustage, x(j,:)] = gf_fft_decode_qec(x(j,:),f(:,j),p,m,omega,lengths(2:end),N);
        
        x_next_stage_order = reshape(x(j,:),lengths(2),prod(lengths(3:end)));
        x(j,:) = reshape(x_next_stage_order',1,[]);
        
        % Multiply with gamma^(i'*k'')
        gamma = mod(omega*(N/(size(x,1)*size(x,2))),p^m-1);  %omega*(N/(len*no_of_blocks))
        for k = 0:(size(x,2)-1)
            i = j-1;
            multiplier = mod(gamma*i*k,p^m-1);
            x(j,k+1) = gfmul(x(j,k+1),multiplier,field);
        end
        
        for k = 1:no_of_blocks
            if (nnz(y(:,k) == 1/p^m) <= j) % #outputs erased <= #inputs known
                inp = x(1:j,k);
                [x(:,k), ~] = forney_fft(y(:,k),inp,p,m,mod(omega*(N/len),p^m-1),field);
            end
        end
        u(:,j) = ustage;   % Accumulate decoded inputs
    end
    u = reshape(u',1,[]);
    u = u';              % Pass on "big code's" inputs to the next stage
    for k = 1:no_of_blocks
        x(:,k) = gf_ifft(x(:,k),p,m,[],length(x(:,k)));  % Return output to the next stage
    end
    if (N == prod(lengths))
        x = reshape(x',1,[]);
        x = x';
    else
        x = reshape(x,1,[]);
    end
end
return
end

function [inp_est, r_est] = forney_fft(r,inp,p,m,omega,field)
% Run Forney's algorithm to correct 'r' based on available data in
% 'inp'; 'inp' contains a few known values that should be obtained
%  when 'r' is evaluated at corresponding initial powers of
%  omega^-1, starting from 0 i.e. @ (omega^-1)^0,(omega^-1)^1 etc.

%disp(omega);
omega = gfdiv(0,omega,field);
erased_indices = find(r==1/p^m);  % Erased indices will have value = 1/q
r_est = r';
% b is the smallest exponent of omega_inv that is a root of g(x)
b = find(inp~=1/p^m,1,'first')-1;

if (~isempty(erased_indices))
    % Compute erasure locator polynomial
    lambda = 0;
    for l = 1:length(erased_indices)
        X_l = gfpow(omega,erased_indices(l)-1,field);  % Xl = omega^i_l
        fact = [0 gfsub(-Inf,X_l,field)];     % fact = (1-X_l*x)
        lambda = gfconv(lambda,fact,field); % lambda(x)=PI((1-X_l*x))
    end
    r_est(erased_indices) = -Inf;         % 0 in exponential format
    %fprintf('Received Word:');
    %disp(r);
    
    no_of_syndromes = length(find(inp~=1/p^m));  % Number of known evaluations
    
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
        
        % Error evaluator polynomial
        xt = [(-Inf)*ones(1,no_of_syndromes) 0];   % x^(2t)
        [~, err_eval] = gfdeconv(gfconv(S,lambda,field),xt,field);
        
        % Find lamba'(x)
        lambda_d = gf_diff_poly(lambda,field);
        
        omega_inv = gfdiv(0,omega,field);
        
        %fprintf('\nErasure Values:');
        % Compute error values and correct r(x)
        for k = 1:length(erased_indices)
            X_k_inv = gfpow(omega_inv,erased_indices(k)-1,field);  % Xl^-1 = (omega^(-i_l))
            err_eval_k = gf_eval_poly(err_eval,X_k_inv,field);
            lambda_d_k = gf_eval_poly(lambda_d,X_k_inv,field);
            err_val_k = gfsub(-Inf,gfdiv(err_eval_k,lambda_d_k,field),field);
            if (b == 0)
                X_k = gfdiv(0,X_k_inv,field);
                err_val_k = gfmul(err_val_k,X_k,field);
            end
            r_est(erased_indices(k)) = gfsub(r_est(erased_indices(k)),err_val_k,field);
        end
        %fprintf('Decoded Codeword:');
        %disp(r_est);
        
        % Compute the input based on corrected r(x)
        inp_est = gf_fft(r_est,p,m,[],length(r_est));
    else
        inp_est = gf_fft(r_est,p,m,[],length(r_est));   % Will return garbage estimate
        %inp_est = (1/p^m)*ones(len_r,1);
    end
else
    inp_est = gf_fft(r_est,p,m,[],length(r_est));   % Exact inputs since no erasures
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