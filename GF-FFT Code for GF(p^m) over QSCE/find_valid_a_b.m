clc
close all
clear all

a_max = 1000;
b_max = 1000;
valid_a_b = zeros(1,2);
index = 0;

for a = 2:a_max
    for b = 2:b_max
        N = a*b;
        violate = 0;
        for i = 0:(N-1)
            s_ba = b*mod(i,a) + floor(i/a);
            s_ab = a*mod(i,b) + floor(i/b);
            if (floor(s_ba/b)*mod(s_ba,b) ~= floor(s_ab/b)*mod(s_ab,b))
                violate = 1;
                break;
            end
        end
        if (violate)
            fprintf('\na = %d, b = %d, Invalid.\n',a,b);
        else
            fprintf('\na = %d, b = %d, NO Violations!\n',a,b);
            index = index + 1;
            valid_a_b(index,:) = [a b];
        end
    end
end

fprintf('\n%d valid pairs (a,b).\n\n',index);