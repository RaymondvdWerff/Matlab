n = 200;
A = rand(n);
A = A + A';

format short
m = 10;
V = zeros(n,m);
T = zeros(m,m);

v = rand(n,1);
v = v/norm(v);
V(:,1) = v;
w_prime = A*V(:,1);
T(1,1) = w_prime'*V(:,1);
w = w_prime - T(1,1)*V(:,1);

for j = 2:m
    T(j-1,j) = norm(w);
    if T(j-1,j) == 0
        disp('error');
    end
    T(j,j-1) = T(j-1,j);
    V(:,j) = w/T(j-1,j);
    w_prime = A*V(:,j);
    T(j,j) = w_prime'*V(:,j);
    w = w_prime - T(j,j)*V(:,j) - T(j-1,j)*V(:,j-1);
    for i = 1:j-1
       w = w - w'*V(:,i)*V(:,i); 
    end
end

Qtot = eye(m);

for i = 1:10
    [Q,R] = qr(T);
    Qtot = Qtot*Q;
    T = R*Q;
end

[vec1,val1] = eigs(A,1,'LM');
[~,index] = max(abs(vec1(:)));
vec1 = vec1*sign(vec1(index));

val1
T(1)

vec2 = V*Qtot(:,1);
[~,index] = max(abs(vec2(:)));
vec2 = vec2*sign(vec2(index));

min(vec1-vec2)
max(vec1-vec2)
%{
m = n;
x = rand(n,1);
q = x/norm(x);
Q = [q];
r = A*q;
a(1) = q'*r;
T(1,1) = a(1);
r = r - a(1)*q;
b(1) = norm(r);

for j = 2:m
    v = q;
    q = r/b(j-1);
    Q = [Q,q];
    r = A*q - b(j-1)*v;
    a(j) = q'*r;
    T(j,j) = a(j);
    r = r-a(j)*q;
    b(j) = norm(r);
    T(j-1,j) = b(j);
    T(j,j-1) = b(j);
    if b(j) == 0
        disp('special case');
        b = b(1:j-1);
        break
    end       
end
Q*Q'
eigs(A,4,'LM')
eigs(T,4,'LM')
%}
