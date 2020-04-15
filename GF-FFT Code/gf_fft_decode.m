function [u,x] = gf_fft_decode(y,f,q,omega,stage,N)
% y = bit APP matrix from channel in output order;
%     each column signifies a particular bit; each row signifies a GF value
% f = input a priori probs matrix in input decoding order;
%     each column signifies a particular bit; each row signifies a GF value
% x = output hard decision probability matrix in output order
% u = input hard decisions in input order
% omega = primitive Nth root of unity in GF(q) where 'q' is prime

len = size(y,2);
if (N == len)
    for i = 2:(q-1)   % Find the smallest primitive Nth root of unity
        % a = mod(i.^[1:N],q); Not suitable because the powers can shoot up
        a(1) = i;
        for j = 2:N
            a(j) = mod(a(j-1)*i,q);
        end
        if (any(a==1) && find(a==1,1,'first') == N)
            omega = i;  % primitive Nth root of unity in GF(q); 'q' prime
            fprintf('%s%d\n','The primitive root was chosen to be ',omega);
            break;
        end
    end
    bin_fft_order = fliplr(dec2bin(0:(N-1),log2(N)));
    fft_order = bin2dec(bin_fft_order) + 1;
    y = y(:,fft_order);
    f = f(:,fft_order);
    stage = zeros(1,log2(N)-1);
end
% Recurse down to length 1
if (len == 1)
    if (all(f == (1/q)*ones(q,1)))
        % If info bit, make hard decision based on observation
        u = find(y == max(y),1,'first') - 1;
        x = zeros(q,1);   x(u+1,1) = 1;  % Change probs. based on hard decision
    else
        % If frozen, use frozen bit
        x = f; u = find(f == 1,1,'first') - 1;
    end
else
    
    % twiddle = mod(omega^(N/2),q);
    twiddle = 1;
    for j = 1:(N/2)
        twiddle = mod(twiddle*omega,q);
    end
    
    mul = mod((q-1)*twiddle,q);  % mul = -omega^(N/2) = (q-1)*omega^(N/2)
    val = gfdiv(1,mod(1+mul,q),q);  % val = (1-omega^(N/2))^(-1) in GF(q)
    % mod(1+mul,q)= 2 always since omega^(N/2)= -1 always
    % omega^N = 1 => omega^(N/2)= +-1. But, since omega is a primitive Nth
    % root of unity, omega^k ~= +1 for k=1,2,...,N-1. Hence omega^(N/2)= -1
    
    % Compute soft mapping back one stage
    u1est = permute_prob(cnop(y(:,2:2:end),permute_prob(y(:,1:2:end),mul,q),q),val,q);
    
    half_len = len/2;          % Same as size(u1est,2)
    stage_index = N/half_len;  % Stage index - CRUCIAL PARAMETER
    
    if (stage_index == 2)
        base_index = 1;
    else
        base_index = stage(log2(stage_index) - 1);
    end
    
    stage(log2(stage_index)) = base_index;               % First time at this stage
    
    % factors = mod(omega.^[0:1:(N/2-1)],q);
    factors(1,1) = 1;
    for j = 2:(N/2)
        factors(1,j) = mod(factors(1,j-1)*omega,q);
    end
    
    current_factors = ones(1,stage_index);
    if (half_len > 1)
        current_factors = [current_factors factors(1,1:(half_len/2):end)];
        current_factors = repmat(current_factors,1,(half_len/2));
    end
    permute_factors = current_factors(1,base_index:stage_index:end);
    permute_factors_inv = gfdiv(ones(1,length(permute_factors)),permute_factors,q);
    
    u1est_permuted = permute_cols(u1est,permute_factors_inv,q);
    
    % R_N^T maps u1est to top polar code
    [uhat1,u1hardprev] = gf_fft_decode(u1est_permuted,f(:,1:(len/2)),q,omega,stage,N);
    u1hardprev = permute_cols(u1hardprev,permute_factors,q);
    
    % Using u1est and x1hard, we can estimate u2
    u2guess1 = cnop(permute_prob(u1hardprev,q-1,q),y(:,1:2:end),q);
    u2guess2 = permute_prob(cnop(y(:,2:2:end),permute_prob(u1hardprev,q-1,q),q),gfdiv(1,twiddle,q),q);
    u2est = vnop(u2guess1,u2guess2,q);
    
    base_index_2 = base_index + (stage_index/2);
    stage(log2(stage_index)) = base_index_2;             % Second time at this stage
    
    permute_factors_2 = current_factors(1,base_index_2:stage_index:end);
    permute_factors_2_inv = gfdiv(ones(1,length(permute_factors_2)),permute_factors_2,q);
    
    u2est_permuted = permute_cols(u2est,permute_factors_2_inv,q);
    
    % R_N^T maps u2est to bottom polar code
    [uhat2,u2hardprev] = gf_fft_decode(u2est_permuted,f(:,(len/2+1):end),q,omega,stage,N);
    u2hardprev = permute_cols(u2hardprev,permute_factors_2,q);
    
    % Tunnel u decisions back up. Compute and interleave x1,x2 hard decisions
    u = [uhat1 uhat2];
    x = reshape([cnop(u1hardprev,u2hardprev,q); cnop(u1hardprev,permute_prob(u2hardprev,twiddle,q),q)],q,[]);
    xhat = zeros(1,size(x,2));
    for i = 1:size(x,2)
        xhat(1,i) = find(uint32(x(:,i)) == 1,1,'first') - 1;
    end
    
    if (N == len)
        u = u(fft_order);
    end
end
return
% Multiply probability vector by a GF(q) value; induces a permutation
% Let X be a R.V. signifying the value at a particular bit position
% Pr(X = a)=p_{a} =>For b in GF(q),Pr(b*X = a)=Pr(X=inv(b)*a)=p_{inv(b)*a}
% Thus, the probability subscripts are permuted in GF(q)
    function pr = permute_prob(p,b,q)
        if (all(b == 1))
            pr = p;
        else
            binv = gfdiv(1,b,q);
            pr = zeros(q,size(p,2));
            for iter = 1:size(p,2)
                p_indices = 0:(q-1);
                new_indices = mod(binv*p_indices,q);
                pr(:,iter) = p(new_indices+1,iter);
            end
        end
        return
    end
% Permute each column with the corresponding permute factor
    function pr = permute_cols(p,perm_fact,q)
        if (all(perm_fact == 1))
            pr = p;
        else
            pr = zeros(q,size(p,2));
            for iter = 1:size(p,2)
                pr(:,iter) = permute_prob(p(:,iter),perm_fact(iter),q);
            end
        end
        return
    end
% Check-node operation in probability domain
    function z = cnop(w1,w2,q)
        z = zeros(q,size(w1,2));
        for iter = 1:size(w1,2)
            z(:,iter) = cconv(w1(:,iter),w2(:,iter),q);
        end
        return
    end
% Bit-node operation in probability domain
    function z = vnop(w1,w2,q)
        z = zeros(q,size(w1,2));
        for iter = 1:size(w1,2)
            z(:,iter) = w1(:,iter).*w2(:,iter);
            z(:,iter) = z(:,iter)./sum(z(:,iter));  % re-normalize
        end
        return
    end
end