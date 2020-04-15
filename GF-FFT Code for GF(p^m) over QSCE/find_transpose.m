function Mt = find_transpose(M)
% M needs to be a square matrix

n = size(M,1);
E = eye(n);
Mt = zeros(n);

N = n^2;
inds = 1 + reshape(reshape(0:(N-1),n,n)',1,[]);
inds2 = sub2ind([N,N],1:N,inds);
S = zeros(N);
S(inds2) = 1;
for i = 1:n
    Mt = Mt + kron(E(:,i)',E) * S * kron(E,M) * E(:) * E(:,i)';
end

end