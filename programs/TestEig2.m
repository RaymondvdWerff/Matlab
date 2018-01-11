q = 2;
X = 40;
tol = 1e-6;
temp = 1.07;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;

Qsq = sqrtm(ones(q)+(exp(1/temp)-1)*eye(q));
A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
[C,T] = beginmatrices(Qsq,A,X);

[C,~,~] = Potts_CTM(A,C,T,q,X,tol,10,temp);
C1 = ncon({C,C},{[-1,1],[1,-2]});
C1 = reshape(C1,X^2,1);
%C1 = rand(X);C1 = C1+C1';C1 = reshape(C1,X^2,1);
tic    
C1 = arnoldi(T,T,C1,@mult1,min(X^2,20),tol);
C1 = reshape(C1,X,X)/max(abs(C1(:)));
toc

T2 = ncon({T,T},{[-3,1,-1],[-4,1,-2]});
T2 = reshape(T2,X^2,X^2);

tic
[C2,~] = eigs(T2,1,'LM');
C2 = reshape(C2,X,X)/max(abs(C2(:)));
toc

tic
[C3,~] = eigs(@(x)mult1(T,T,x),X^2,1,'LM');
C3 = reshape(C3,X,X)/max(abs(C3(:)));
toc

%delta12 = ncon({abs(C1-C2),abs(C1-C2)},{[1,2],[1,2]})
%delta13 = ncon({abs(C1-C3),abs(C1-C3)},{[1,2],[1,2]})
%delta23 = ncon({abs(C2-C3),abs(C2-C3)},{[1,2],[1,2]})

[U,s2,~] = svd(C1);
C = U*sqrt(s2)*U';
C = C./max(abs(C(:)));

CT = ncon({C,T},{[-1,1],[1,-2,-3]});
CT = reshape(CT,q*X,X);

[U,s,V] = svd(CT,0);
Tl = U*V';
Tl = reshape(Tl,[X,q,X]);
C_prime = V*s*V';
C_prime = C_prime./max(abs(C_prime(:)));
    
delta = sqrt(ncon({C-C_prime,conj(C-C_prime)},{[1,2],[1,2]}))
    
function v = arnoldi(M1,M2,b,mult,m,tol)
    
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

function y = mult1(M1,M2,C)
    si = size(M1);
    C = reshape(C,si(1),si(1));
    y = ncon({M1,M2,C},{[1,3,-2],[2,3,-1],[2,1]});
    y = reshape(y,si(1)^2,1);
end

function [C0,T0] = beginmatrices(Qsq,A,X)
    q = size(A,1);
    spin1_2D = zeros(q,q);spin1_2D(1,1) = 1;
    spin1_3D = zeros(q,q,q);spin1_3D(1,1,1)=1;
    C0 = ncon({spin1_2D,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T0 = ncon({spin1_3D,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    
    while size(T0,1) < X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TA = ncon({T0,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TA},{[-1,1,2],[1,2,-2,-3,-4]});
        C0 = reshape(M,q*size(T0,1),q*size(T0,1));
        T0 = reshape(TA,[q*size(T0,1),q,q*size(T0,1)]);
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
    
    if size(T0,1) > X
        [U,~,~] = svd(C0);
        U_til = U(:,1:X);
        C0 = ncon({C0,U_til,U_til},{[1,2],[1,-1],[2,-2]});
        T0 = ncon({T0,U_til,U_til},{[1,-2,2],[1,-1],[2,-3]});
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
end