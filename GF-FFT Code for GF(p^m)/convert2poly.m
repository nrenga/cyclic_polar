function [ a_x ] = convert2poly( a )
% convert2poly: To convert the given coff. into a polynomial string
%               for display purposes

if (a(1,1) == 0)
    a_x = '';
else
    a_x = num2str(a(1,1));
end

for i=2:length(a)
    if (a(1,i) ~= 0)
        if (a(1,i) ~= 1)
            append = num2str(a(1,i));
        else
            append = '';
        end
        if (i ~= 2)
            append = strcat(append,'x^',num2str((i-1)));
        else
            append = strcat(append,'x');
        end
        if (strcmp(a_x,'') == 1)
            a_x = strcat(a_x,append);
        else
            a_x = strcat(a_x,' +',' ',append);
        end
    end
end

if (strcmp(a_x,'') == 1)
    a_x = '0';
end
end