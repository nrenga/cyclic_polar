function [ f, I ] = gf_fft_design_qsce(lengths,perror,perasure,d,p,m)
% lengths is a vector of factors of the blocklength N
% Design FFT code of length prod(lengths) for BEC(e) and target block error rate d
% Over GF(p^m) where p is prime and N|(p^m-1)
% Generate virtual channel erasure probabilities
% f = input a priori probs matrix in input order;
%     each column signifies a particular bit; each row signifies a GF value

q = p^m;
N = prod(lengths);  % The actual blocklength of the code
M = 10^4;   % # iterations for Monte Carlo DE

%filenm = strcat('N_',num2str(N),'_GF_',num2str(q),'_QSCE_',num2str(perror),'_',num2str(perasure),'.mat');
%load(filenm);

inpError = perror*ones(1,N);  % Assumed to be in output order
inpErasure = perasure*ones(1,N);  % Assumed to be in output order

% Assume N = n'*n'' is given as N = [n' n'']; n'' could be split into
% multiple factors
fft_order = [0:(N-1)]';
for i = 1:length(lengths)   % Index for the small-block length in the current stage
    remaining_len = prod(lengths((i+1):end));
    no_of_blocks = prod([lengths(1:(i-1)) remaining_len]);
    len = lengths(i);
    inpError = reshape(inpError',1,[]);
    inpError = reshape(inpError,len,no_of_blocks);  % Each column corresponds to the error probs. of a block
    inpErasure = reshape(inpErasure',1,[]);
    inpErasure = reshape(inpErasure,len,no_of_blocks);  % Each column corresponds to the erasure probs. of a block
    [~,indices1,~] = unique(inpError(1,:));  % same as unique(inpErasure(1,:))
    [~,indices2,~] = unique(inpErasure(1,:));
    indices = union(indices1,indices2);
    for ind = 1:length(indices)
        j = indices(ind);
        fprintf('\nlen = %d, perror = %1.4f, perasure = %1.4f\n',len,inpError(1,j),inpErasure(1,j));
        % The next line and the pair of subsequent lines are mutually
        % exclusive. If the file, filenm, exists then comment next line
        [inpError(:,j), inpErasure(:,j)] = monte_carlo_de2(len,p,m,inpError(1,j),inpErasure(1,j),M); % The input erasure probs. for each block
        inpError(:,j:(j+remaining_len-1)) = repmat(inpError(:,j),1,remaining_len);
        inpErasure(:,j:(j+remaining_len-1)) = repmat(inpErasure(:,j),1,remaining_len);
        % Comment next two lines if the file, filenm, does not exist
        %inpError(:,j) = E{i,1}(:,j);
        %inpErasure(:,j) = E{i,2}(:,j);
    end
    fprintf('\ninpError = \n');
    disp(inpError);
    fprintf('\ninpErasure = \n');
    disp(inpErasure);
    
    % Ignoring reordering of E since the probabilities will not change even
    % after that. Reordering here will be painful (but unnecessary)
    tempError = inpError;
    inpError = [];
    tempErasure = inpErasure;
    inpErasure = [];
    for j = 1:(no_of_blocks/remaining_len)
        inpError = [inpError; tempError(:,((j-1)*remaining_len +1):(j*remaining_len))];
        inpErasure = [inpErasure; tempErasure(:,((j-1)*remaining_len +1):(j*remaining_len))];
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
        disp(fft_order);
    end
    %disp(E);
end
C = find_capacity_qsce(q,inpError,inpErasure);
C(isnan(C)) = 0;
%fprintf('The input capacities in order of indices to blocks are: \n');
%disp(C');
fft_order = reshape(fft_order,1,[]);
[~, actual_loc] = sort(fft_order);
C = C(actual_loc);
inpError = inpError(actual_loc);
inpErasure = inpErasure(actual_loc);
fprintf('The input capacities in ascending order of input indices are: \n');
disp(C');
fprintf('The input erasure probs. in ascending order of input indices are: \n');
disp(inpErasure');

avg_output_capacity = find_capacity_qsce(q,perror,perasure);
avg_input_capacity = mean(C);
fprintf('\nAverage output channel capacity = %f\n',avg_output_capacity);
fprintf('\nAverage input channel capacity = %f\n',avg_input_capacity);

inpCorrupt = inpError + inpErasure;
[asc_err, ord] = sort(inpCorrupt);
CE = cumsum(asc_err);
k = nnz(CE <= d);    % No. of info bits
fprintf('\nRate = %d/%d = %f\n',k,N,k/N);
f = (-Inf)*ones(N,1);
I = ord(1:k);
f(I) = 1/q;
fprintf('The positions of information are: \n');
disp(I'-1);

% filenm = strcat('N_',num2str(N),'_GF_',num2str(q),'_QSCE_',num2str(perror),'_',num2str(perasure),'.mat');
% save(filenm);
end