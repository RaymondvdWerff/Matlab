function [T,C,iter] = Potts_CTM(T,C,A,A1,q,X,tol,temp)
    
    delta = tol + 1;    
    iter = 0;
    maxiter = 500;
    sv = randn(X);
    
    while delta > tol
        
        if iter > 0
            A = A1;
        end
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,sv_new,~] = tensorsvd(M,[1,2],[3,4],X);
        sv_new = sv_new/max(sv_new(:));
    
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/2;
        C = C/max(C(:));
        
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/2;
        T = T/max(T(:));
        
        delta = 0;
        for k = 1:sqrt(numel(sv_new))
            delta = delta + abs(sv(k,k)-sv_new(k,k));
        end
        
        sv = sv_new;
        iter = iter + 1;
        if iter > maxiter
            disp(['Potts_CTM not converged at T = ' num2str(temp)]);
            break
        end            
    end       
end
 