function [Al,C] = LeftOrthonormalize4(A,C,tol,temp)
    
    iter = 0;
    maxiter = 500;
    X = size(A,1);
    
    [~,C] = qrpos(C);
    C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

    Ac = ncon({C,A},{[-1,1],[1,-2,-3]});
    Cold = C;

    si = size(Ac);
    Ac = reshape(Ac,prod(si(1:2)),si(3));
    [Al,C] = qrpos(Ac);    
    Al = reshape(Al,[si(1),si(2),si(3)]);
    
    lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
    C = C./lambda;
    
    delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));

    while delta > tol
        
        C = arnoldi(A,Al,reshape(C,X^2,1),min(X^2,15),tol);
        C = reshape(C,X,X);
        
        [~,C] = qrpos(C);
        C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

        Ac = ncon({C,A},{[-1,1],[1,-2,-3]});
        Cold = C;

        si = size(Ac);
        Ac = reshape(Ac,prod(si(1:2)),si(3));
        [Al,C] = qrpos(Ac);
        
        Al = reshape(Al,[si(1),si(2),si(3)]);

        lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
        C = C./lambda;
        
        delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));
        
        iter = iter +1;
        if iter > maxiter
            disp(['LeftOrthonormalize4 not converged at T = ' num2str(temp)]);
            break;
        end
    end
end

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