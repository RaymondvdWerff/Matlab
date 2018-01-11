X = 3;
q = 2;
tol = 10^(-4);

T0 = rand(X,q,X);
T0 = T0 + permute(T0,[3,2,1]);

T2 = ncon({T0,conj(T0)},{[-3,1,-1],[-4,1,-2]});
T2 = reshape(T2,X^2,X^2);

[C2,~] = eigs(T2,1,'LM');
C2 = reshape(C2,X,X);

[U,s2,~] = svd(C2);
C = U*sqrt(s2)*U';
C = C./max(abs(C(:)));

CT = ncon({C,T0},{[-1,1],[1,-2,-3]});
CT = reshape(CT,q*X,X);

[U,s,V] = svd(CT,0);
Tl = U*V';
Tl = reshape(Tl,[X,q,X]);
C1 = V*s*V';
delta = sqrt(ncon({C-C1,conj(C-C1)},{[1,2],[1,2]}));

for i = 1:10
    
    C1 = arnoldi(T0,Tl,reshape(C1,X^2,1),min(X^2,15),tol);
    C1 = reshape(C1,X,X);
    
    [U,s,V] = svd(C1);
    Q = U*V';
    C_prime = V*s*V';
    
    CT = ncon({C_prime,T0},{[-1,1],[1,-2,-3]});
    CT = reshape(CT,q*X,X);

    [U,s,V] = svd(CT,0);
    Tl = U*V';
    Tl = reshape(Tl,[X,q,X]);
    C1 = V*s*V';
    C1 = C1./max(abs(C1(:)));
    delta = sqrt(ncon({C-C1,conj(C-C1)},{[1,2],[1,2]}));
    
end

%C1./max(abs(C1(:)))-C./max(abs(C(:)))


M1 = ncon({Tl,C},{[-1,-2,1],[1,-3]});
M2 = ncon({C,T0},{[-1,1],[1,-2,-3]});
M1 = M1./max(abs(M1(:)));
M2 = M2./max(abs(M2(:)));
M1-M2


function v = arnoldi(M1,M2,b,m,tol)
    
    H = zeros(m+1,m);
    V = zeros(size(b,1),m+1);
    V(:,1) = b/norm(b);
    
    for j = 1:m
        w = mult(M1,M2,V(:,j));
        for i = 1:j
            H(i,j) = w'*V(:,i);
            w = w - H(i,j)*V(:,i);
        end
        
        H(j+1,j) = norm(w);
        
        if H(j+1,j) < tol
            H = H(1:j,1:j);
            break
        end
        
        V(:,j+1) = w./H(j+1,j);
    end
    
    H = H(1:j,1:j);
    V = V(:,1:j);
    [vector,~] = eigs(H,1,'LM');
    v = V*vector;
    
end

function y = mult(M1,M2,C)
    si = size(M1);
    C = reshape(C,si(1),si(1));
    y = ncon({conj(M2),C,M1},{[1,3,-1],[1,2],[2,3,-2]});
    y = reshape(y,si(1)^2,1);
end