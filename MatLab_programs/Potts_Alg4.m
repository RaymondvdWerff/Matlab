function [T,C,iter,tictoc] = Potts_Alg4(T,C_notused,A,q,X,tol,temp)
    
    delta = tol + 1;
    iter = 0;
    maxiter = 1000;
    opts.isreal = false;
    sold = rand(X);
    
    tic;
    while delta > tol
        
        [Tl,C] = LeftOrthonormalize5(T,1e-6,temp);
        [~,s,~] = svd(C);s = s./max(s(:));
        
        T = arnoldi(permute(Tl,[3,2,1]),A,reshape(T,q*X^2,1),min(q*X^2,15),tol);
        T = reshape(T,[X,q,X]);
        T = T + permute(T,[3,2,1]);
        
        delta = sum(sum(abs(s-sold))); 
        
        iter = iter + 1;
        if iter > maxiter
            disp(['Potts_Alg4 not converged at T = ' num2str(temp)]);
            break
        end
        sold = s;
    end
    tictoc = toc;
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