function [partial_order, violations, permutations, inds_i, I, I_coords, O, O_coords, E, Ecomp, bases] = channel_partial_orders(lengths,e)
% 'lengths' should contain FFT lengths starting from the o/p stage

% clc
% close all
% clear all
% 
% lengths = [2 5 3];
% q = 31;
% 
% [f, I, O, E] = gf_fft_design_qec(lengths,0.5,10^7,q);
[ I, E, Ecomp, O] = gf_fft_qec_precise_order(lengths,e);
I_coords = gf_fft_expand(lengths,I, 0);
O_coords = gf_fft_expand(lengths,O, 0);

violations = cell(1,6);
v_ind = 0;
n = length(lengths);
N = prod(lengths);

bases = zeros(1,n);
for i = 1:n
    bases(i) = prod(lengths(1:(i-1)));
end
bases = sort(bases,'descend');

inds = 0:(N-1);
partial_order = cell(1,1);
order_index = 0;
no_of_perm = 0;

fprintf('\nPartial Orders:\n');
for i = 2:(N-1)
    if ((gcd(i,N) == 1))% && (mod(i^2,N) == 1))
        %fprintf('\nCase i = %d:\n',i);
        no_of_perm = no_of_perm + 1;
        permutations(no_of_perm) = i;
        inds_i{no_of_perm,1} = mod(i*inds,N);
        O_permuted_inds = inds_i{no_of_perm,1}(O+1);
        accounted_inds = [];
        for j = N:-1:2
            if (isempty(find(accounted_inds == O(j),1)))
                accounted_inds = [accounted_inds O(j)];
                perm_index = find(O==O_permuted_inds(j),1);
                if (perm_index < j)
                    if (all( ismember(O(1:perm_index-1),inds_i{no_of_perm}(O(1:j-1)+1)) ))
                        order_index = order_index + 1;
                        partial_order{order_index,1} = O(j);
                        accounted_inds = [accounted_inds O_permuted_inds(j)];
                        a = find(I==O(j));
                        b = find(I==O_permuted_inds(j));
                        %fprintf('\n%s\n',mat2str(partial_order{order_index}));
                        if (b < a)
                            v_ind = v_ind + 1;
                            violations{1}(v_ind) = I(a);
                            violations{2}(v_ind,1:n) = I_coords(a,:);
                            violations{3}(v_ind,:) = [E(I(a)+1) Ecomp(I(a)+1)];
                            violations{4}(v_ind) = I(b);
                            violations{5}(v_ind,1:n) = I_coords(b,:);
                            violations{6}(v_ind,:) = [E(I(b)+1) Ecomp(I(b)+1)];
                        end
                    end
                end
            end
        end
    end
end
%fprintf('\n');

fprintf('\n# of pair-orderings = %d and there were %d violations.\n\n',order_index,v_ind);

end