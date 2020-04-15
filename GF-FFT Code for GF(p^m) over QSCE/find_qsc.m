function [perror] = find_qsc(rate,q)
% This function computes QSC(perror) such that rate <= capacity of
% QSC(perror)

precision = 0.001;
for p = 1:(-precision):precision
    H = p*log_base(q-1,q) - p*log_base(p,q) - (1-p)*log_base(1-p,q);
    C = 1-H;
    if (C >= rate)
        perror = p;
        break;
    end
end
end