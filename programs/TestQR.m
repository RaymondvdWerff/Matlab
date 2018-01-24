D = 10;
A = rand(D);
A = A+A';

format long 

[V,~] = eigs(A,1,'LM');

[~,index] = max(abs(V(:)));
V = V*sign(V(index));

Qtot = eye(D);

for i = 1:14
    [Q,R] = qr(A);
    Qtot = Qtot*Q;
    A = R*Q;
end
vec1 = Qtot(:,1);
[~,index] = max(abs(vec1(:)));
vec1 = vec1*sign(vec1(index));

min(V-vec1)
max(V-vec1)
