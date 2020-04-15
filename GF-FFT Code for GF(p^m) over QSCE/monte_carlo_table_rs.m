function [ inp_channel_table ] = monte_carlo_table_rs( len,p,m,M )
% Perform Monte Carlo Simulations to determine DE results for the given
% code on the channel QSCE(perror,perasure) using a hard decoder

q = p^m;
verbose = 0;

% Row index - # errors, Col index - # erasures
% Entry at (i,j) = [<inperror vector> <inperasure vector>]
inp_channel_table = cell(len+1, len+1);

error_erasure_index = zeros(nchoosek(len+2,2),2);
ind = 1;
for no_of_errors = 0:len
    for no_of_erasures = 0:(len-no_of_errors)
        error_erasure_index(ind,:) = [no_of_errors no_of_erasures];
        ind = ind + 1;
    end
end

for ind = 1:nchoosek(len+2,2)
    
    no_of_errors = error_erasure_index(ind,1);
    no_of_erasures = error_erasure_index(ind,2);
    
    inperror = zeros(len,1);
    inperasure = zeros(len,1);
    
    u = randi(q,1,len)-2;
    u(u == -1) = -Inf;
    
    v = gf_ifft(u,p,m,[],len);
    
    %fprintf('\nlen = %d\n',len);
    for i = 1:M
        clc
        fprintf('\nindex = %d of %d\n',ind,nchoosek(len+2,2));
        fprintf('\n# errors = %d, # erasures = %d\n',no_of_errors,no_of_erasures);
        fprintf('\nlen = %d, i = %d\n',len,i);
        
        r = v;
        
        % QSCE Channel with given #errors, #erasures
        % Doesn't matter which indices are errored or erased
        error_indices = 1:no_of_errors;
        erased_indices = (no_of_errors+1):(no_of_errors + no_of_erasures);
        
        symbols = [-Inf 0:(q-2)];
        for j = 1:no_of_errors
            new_symbols = symbols;
            new_symbols(new_symbols == v(error_indices(j))) = [];
            r(error_indices(j)) = new_symbols(randi(q-1));
        end
        
        r(erased_indices) = 1/q;
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
    inp_channel_table{no_of_errors+1, no_of_erasures+1} = [inperror inperasure];
    filenm = strcat('len_',num2str(len),'_M_',num2str(M),'.mat');
    save(filenm,'inp_channel_table');
end
disp(error_erasure_index);

filenm = strcat('len_',num2str(len),'_M_',num2str(M),'.mat');
save(filenm,'inp_channel_table');
end