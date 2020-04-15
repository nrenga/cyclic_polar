function [violations, swaps, I, I_coords, O, O_coords, E, Ecomp, bases] = check_bitsig_order(lengths,e)
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

violations = cell(1,4);
swaps = cell(1,4);
n = length(lengths);

bases = zeros(1,n);
for i = 1:n
    bases(i) = prod(lengths(1:(i-1)));
end
bases = sort(bases,'descend');

lengths = fliplr(lengths);
v_ind = 0;
combos = nchoosek(1:n,2);
no_of_swaps = 0;
for i = 1:length(I)
    for j = 1:size(combos,1)
        p = combos(j,1);
        q = combos(j,2);
        swap_coord = I_coords(i,:);
        if ((lengths(p) == lengths(q)) && (I_coords(i,p) < I_coords(i,q))) % && (I_coords(i,p) < lengths(q)) && (I_coords(i,q) < lengths(p)))
            swap_coord(p) = I_coords(i,q);
            swap_coord(q) = I_coords(i,p);
            no_of_swaps = no_of_swaps + 1;
            swap_index = swap_coord*bases';
            if (find(I==swap_index) < i)
                v_ind = v_ind + 1;
                violations{1}(v_ind) = I(i);
                violations{2}(v_ind,1:n) = I_coords(i,:);
                violations{3}(v_ind,:) = [E(I(i)+1) Ecomp(I(i)+1)];
                violations{4}(v_ind) = swap_index;      
                violations{5}(v_ind,1:n) = swap_coord;
                violations{6}(v_ind,:) = [E(swap_index+1) Ecomp(swap_index+1)];
            end
            swaps{1}(no_of_swaps) = I(i);
            swaps{2}(no_of_swaps,1:n) = I_coords(i,:);
            swaps{3}(no_of_swaps,:) = [E(I(i)+1) Ecomp(I(i)+1)];
            swaps{4}(no_of_swaps) = swap_index;
            swaps{5}(no_of_swaps,1:n) = swap_coord;
            swaps{6}(no_of_swaps,:) = [E(swap_index+1) Ecomp(swap_index+1)];
        end
    end
end
fprintf('\n# of swaps = %d and there were %d violations.\n\n',no_of_swaps,v_ind);

end