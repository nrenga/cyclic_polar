function [ inperror, inperasure ] = monte_carlo_de( len,p,m,perror,perasure,M )
% Perform Monte Carlo Simulations to determine DE results for the given
% code on the channel QSCE(perror,perasure) using a hard decoder

q = p^m;
verbose = 0;

inperror = zeros(len,1);
inperasure = zeros(len,1);

u = randi(q,1,len)-2;
u(u == -1) = -Inf;

v = gf_ifft(u,p,m,[],len);

%fprintf('\nlen = %d\n',len);
for i = 1:M
    clc
    fprintf('\nlen = %d, i = %d\n',len,i);
    
    r = v;
    
    % QSCE Channel
    pcorrect = 1 - perror - perasure;
    symbols = [1/q -Inf 0:(q-2)];
    for j = 1:len
        pvec = [perasure (perror/(q-1))*ones(1,q)];
        pvec(symbols == v(j)) = pcorrect;
        plen = length(pvec);     % plen = q+1
        pvec = reshape(pvec,plen,1);   % pvec = pvec'
        index = sum(repmat(rand(1),plen,1) > cumsum(pvec),1) + 1;
        r(j) = symbols(index);
    end
    %disp(r);
    
    inp = (1/q)*ones(len,1);
    for j = 1:len
        [inp_est,~] = bm_forney_fft(r,inp,p,m,verbose);
        %disp(r_est);
        if (inp_est(j) == 1/q)
            inperasure(j) = inperasure(j) + 1;
        elseif (inp_est(j) ~= u(j))
            inperror(j) = inperror(j) + 1;
        end
        inp(j) = u(j);
    end
end
inperasure = inperasure/M;
inperror = inperror/M;

end