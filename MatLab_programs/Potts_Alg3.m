function [T1,C,iter] = Potts_Alg3(T1,C_prime,A,A1,q,X,tol,temp)
    
    delta = tol + 1;
    iter = 0;
    maxiter = 1000;
    opts.isreal = false;
    T2 = T1;
    
    while delta > tol
        
        if iter > 0
            A = A1;
        end
        
        [T1l,C] = LeftOrthonormalize4(T1,C_prime,tol,temp);
        T1l = real(T1l);C = real(C);
        
        T2 = arnoldi(permute(T1l,[3,2,1]),A,reshape(T2,q*X^2,1),min(q*X^2,15),tol);
        T2 = reshape(T2,[X,q,X]);
        T2 = T2 + permute(T2,[3,2,1]);
        
        [T2r,C_prime] = LeftOrthonormalize4(permute(T2,[3,2,1]),C,tol,temp);
        T2r = real(T2r);C_prime = real(C_prime);
        T2r = permute(T2r,[3,2,1]);
        C_prime = permute(C_prime,[2,1]);
        
        T1 = arnoldi(T2r,A,reshape(T1,q*X^2,1),min(q*X^2,15),tol);
        T1 = reshape(T1,[X,q,X]);
        T1 = T1 + permute(T1,[3,2,1]);
        
        delta = sqrt(ncon({C-C_prime,conj(C-C_prime)},{[1,2],[1,2]})); 
        
        iter = iter + 1;
        if iter > maxiter
            disp(['Potts_Alg not converged at T = ' num2str(temp)]);
            break
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


function y = mult(M1,M2,x)
    si = size(M1);
    x = reshape(x,[si(1),si(2),si(3)]);
    y = ncon({x,M1,M2,M1},{[2,1,4],[-3,3,2],[3,5,-2,1],[-1,5,4]});
    y = reshape(y,prod(si),1);
end