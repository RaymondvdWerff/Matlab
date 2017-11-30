function [Al,C] = LeftOrthonormalize5(A,tol,temp)
    
    iter = 0;
    miniter = 2;
    maxiter = 500;
    X = size(A,1);
    q = size(A,2);

    A2 = ncon({A,conj(A)},{[-3,1,-1],[-4,1,-2]});
    A2 = reshape(A2,X^2,X^2);

    [C2,~] = eigs(A2,1,'LM');
    C2 = reshape(C2,X,X);

    [U,s2,~] = svd(C2);
    C = U*sqrt(s2)*U';
    C = C./max(abs(C(:)));

    CA = ncon({C,A},{[-1,1],[1,-2,-3]});
    CA = reshape(CA,q*X,X);

    [U,s,V] = svd(CA,0);
    Al = U*V';
    Al = reshape(Al,[X,q,X]);
    C1 = V*s*V';
    
    delta_old = 1;
    for iter = 1:maxiter

        C1 = arnoldi(A,Al,reshape(C1,X^2,1),min(X^2,15),tol);
        C1 = reshape(C1,X,X);

        [~,s,V] = svd(C1);
        C_prime = V*s*V';

        CA = ncon({C_prime,A},{[-1,1],[1,-2,-3]});
        CA = reshape(CA,q*X,X);

        [U,s,V] = svd(CA,0);
        Al = U*V';
        Al = reshape(Al,[X,q,X]);
        C1 = V*s*V';
        C1 = C1./max(abs(C1(:)));
        delta = sqrt(ncon({C-C1,conj(C-C1)},{[1,2],[1,2]}));
              
        if (iter > miniter) && (delta < tol)
            break;
        end
        
        if (delta-delta_old < tol) && (delta < 1e-5)
            break;
        end       
        delta_old = delta;
    end
    if iter == maxiter
            disp(['LeftOrthonormalize5 not converged at T = ' num2str(temp)]);
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