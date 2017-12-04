function [iters,ticks,m] = converge_m_Alg4(T,A0,A1,q,X,temp,maxiter,tol)

    opts.isreal = false;
    tic
    
    for iter = 1:maxiter
        
        if iter == 1
            A = A0;
        end
        if iter > 1
            A = A1;
        end
        
        [Tl,C] = LeftOrthonormalize5(T,1e-6,temp);
        Tl = real(Tl);C = real(C);
        
        T = arnoldi(permute(Tl,[3,2,1]),A,reshape(T,q*X^2,1),min(q*X^2,15),1e-6);
        T = reshape(T,[X,q,X]);
        T = T + permute(T,[3,2,1]);
        
        n = collapse(C,T,A0)/collapse(C,T,A1);
        m(iter) = (q*n-1)/(q-1);
        ticks(iter) = toc;
        
        if iter > 1
            if abs(m(iter)-m(iter-1)) < tol
                break;
            end
        end
    end
    iters = 1:iter;
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