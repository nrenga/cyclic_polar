function coords = gf_fft_expand(lengths, indices, verbose)
% Get coordinates based on the FFT blocks and their ordering
% 'lengths' is a vector of various FFT lengths beginning from the o/p stage
% 'indices' is a vector of all indices for which expansion is desired

l = length(lengths);
bases = zeros(1,l);
for i = 1:l
    bases(i) = prod(lengths(1:(i-1)));
end
bases = sort(bases,'descend');

coords = zeros(length(indices),l);
if (verbose)
    fprintf('\nIndex ---> Coordinates for Lengths = %s, hence Bases = %s: \n',mat2str(fliplr(lengths)),mat2str(bases));
end
for i = 1:length(indices)
    val = indices(i);
    for j = 1:l
        coords(i,j) = floor(val/bases(j));
        val = mod(val, bases(j));
    end
    % Display
    if (verbose && (length(num2str(indices(i))) == 1))
        fprintf(' ');
    end
    if (verbose)
        fprintf('%d ---> %s\n',indices(i),mat2str(coords(i,:)));
    end
end

end