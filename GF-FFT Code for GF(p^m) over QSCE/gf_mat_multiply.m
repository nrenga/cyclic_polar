function C = gf_mat_multiply(A,B,field)

[mA, nA] = size(A);
[mB, nB] = size(B);
C = (-Inf)*ones(mA,nB);  % exponential notation

for i = 1:mA
    for j = 1:nB
        for k = 1:nA
            C(i,j) = gfadd(C(i,j), gfmul(A(i,k),B(k,j),field), field);
        end
    end
end

end