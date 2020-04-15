function [ inperror, inperasure ] = monte_carlo_de2( len,p,m,perror,perasure,M )
    
q = p^m;
%filenm = strcat('len_',num2str(len),'_GF_',num2str(q),'_M_',num2str(M),'.mat');
filenm = strcat('len_',num2str(len),'_M_',num2str(M),'.mat');
load(filenm,'inp_channel_table');

v = (-Inf)*ones(len,1);  % Use all-zero codeword wolog
r = v;

inperror = zeros(len,1);
inperasure = zeros(len,1);

M = 10^5;
for i = 1:M
    no_of_errors = 0;
    no_of_erasures = 0;
    
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
        if (r(j) == 1/q)
            no_of_erasures = no_of_erasures + 1;
        elseif (r(j) ~= v(j))
            no_of_errors = no_of_errors + 1;
        end
    end
    
    inperror = inperror + inp_channel_table{no_of_errors+1,no_of_erasures+1}(:,1);
    inperasure = inperasure + inp_channel_table{no_of_errors+1,no_of_erasures+1}(:,2);
end
inperasure = inperasure/M;
inperror = inperror/M;

end