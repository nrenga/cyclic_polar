function [ f ] = gf_fft_design_qec(lengths,e,d,q)
% N is a vector of co-prime lengths
% Design FFT code of length prod(N) for BEC(e) and target block error rate d
% Over GF(q^m) where q is prime and N|(q^m-1)
% Generate virtual channel erasure probabilities
% f = input a priori probs matrix in input order;
%     each column signifies a particular bit; each row signifies a GF value

N = prod(lengths);  % The actual blocklength of the code
E = e*ones(1,N);  % Assumed to be in output order

% Assume N = n'*n'' is given as N = [n' n'']; n'' could be split into
% multiple factors
fft_order = [0:(N-1)]';
for i = 1:length(lengths)   % Index for the small-block length in the current stage
    remaining_len = prod(lengths((i+1):end));
    no_of_blocks = prod([lengths(1:(i-1)) remaining_len]);
    len = lengths(i);
    E = reshape(E',1,[]);
    E = reshape(E,len,no_of_blocks);  % Each column corresponds to the erasure probs. of a block
    for j = 1:no_of_blocks
        e_prob = E(1,j);         % All values in a column will be equal
        pd = makedist('Binomial','N',len,'p',1-e_prob);  % Success = Not Erased
        E(:,j) = cdf(pd,[(len-1):-1:0]'); % The input erasure probs. for each block
    end    
    % Ignoring reordering of E since the probabilities will not change even
    % after that. Reordering here will be painful (but unnecessary)
    tempE = E;
    E = [];
    for j = 1:(no_of_blocks/remaining_len)
        E = [E; tempE(:,((j-1)*remaining_len +1):(j*remaining_len))];
    end
    % Reorder indices
    if (i < length(lengths))
        old_fft_order = fft_order;
        fft_order = [];
        for j = 1:size(old_fft_order,2)
            for offset = 1:len
                fft_order = [fft_order old_fft_order(offset:len:end,j)];
            end
        end
        %disp(fft_order);
    end
    %disp(E);
end 
fft_order = reshape(fft_order,1,[]);
[asc, actual_loc] = sort(fft_order);
E = E(actual_loc);
fprintf('The input erasure probabilities are: \n');
disp(E');

[asc_err, ord] = sort(E);
CE = cumsum(asc_err);
k = nnz(CE <= d);    % No. of info bits
%fprintf('Rate = %d/%d = %f\n',k,N,k/N);
f = (-Inf)*ones(N,1);
f(ord(1:k)) = 1/q;
fprintf('The positions of information are: \n');
disp(ord(1:k)'-1);
end