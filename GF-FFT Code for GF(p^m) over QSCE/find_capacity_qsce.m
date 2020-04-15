function [C] = find_capacity_qsce(q,perror,perasure)
% This function computes capacity of qSCE(perror,perasure)

C = zeros(length(perror),1);
for i = 1:length(perror)
    pcorrect = 1 - perror(i) - perasure(i);
    H = 0;
    if (perror(i) ~= 0)
        H = perror(i)*log_base(q-1,q) + perror(i)*log_base(pcorrect/perror(i),q) - (1-perasure(i))*log_base(pcorrect/(1-perasure(i)),q);
    end
    C(i) = (1 - perasure(i)) - H;
end
end