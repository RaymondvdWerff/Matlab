%{
q = 2;
X = 4;
tol = 10^(-8);

delta = zeros(q,q,q,q);
for i = 1:q
    delta(i,i,i,i) = 1;    
end

spin = zeros(q,q,q,q);
spin(1,1,1,1) = 1;
spin(2,2,2,2) = -1;

Temp = 1;

Q = [exp(1/Temp) exp(-1/Temp); exp(-1/Temp) exp(1/Temp)];
Qsq = sqrtm(Q);

T = zeros(q,q,q);
T(1,1,1) = 1;
T = ncon({T,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});

a = ncon({delta,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
b = ncon({spin,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});

maxiter = 5000;
minsteps = 4;
sold = zeros(q,1);

for i = 1:maxiter

    N = ncon({T,a},{[-1,-3,1],[1,-2,-4,-5]});
    si_N = size(N);
    N_SVD = reshape(N,prod(si_N(1:2)),prod(si_N(3:5)));
    D = min(X,size(N_SVD,1));

    [U,s,~] = svd(N_SVD,'econ');

    U = U(:,1:D);
    U = reshape(U,[si_N(1),si_N(2),D]);
    s = diag(s);
    s = s(1:D);
    snew = s/max(s);

    T = ncon({conj(U),N,U},{[1,2,-1],[1,2,3,4,-3],[3,4,-2]}); 
    T = T+permute(T,[2 1 3]);
    T = T./max(abs(T(:)));

    if numel(sold) == numel(snew)
        diffs = norm(snew-sold);
    else
        diffs = inf;
    end

    if (diffs<tol)  && (i>minsteps) 
        break;
    end
    sold = snew;

end

M = ncon({T,T},{[-1,-3,1],[-2,-4,1]});
si = size(M);
M = reshape(M,prod(si(1:2)),prod(si(3:4)));
M = (M+M')/2;
%}

M = magic(16);

%[v1,D] = eigs(@(x)(M*x),size(M,1),1,'LM');
%v1
%{
[v1,D,v3] = eig(M);
v3(:,1)

v2 = rand(1,size(M,1));
for j = 1:100
    v2 = v2*M;
    v2 = v2./max(abs(v2(:)));
end
v2./norm(v2)
%}

v = zeros(16,1);
v(1) = 1;

%[V,H] = Arnoldi(M,v,size(M,1),10^(-10));
%[vec1,val1] = eigs(H,1,'LM');
%V*vec1

[vec2,val2] = eigs(M,1,'LM');

vec1 = LeadingEigenvector(M,v,10^(-10))
vec2

