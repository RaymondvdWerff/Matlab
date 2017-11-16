function [T,C,iter] = Potts_CTM2(T,C,A,A1,q,X,tol,temp)
    
    delta = tol + 1;    
    iter = 0;
    maxiter = 5000;
    Cold = rand(X,X);
    
    while delta > tol
        
        if iter > 0
            A = A1;
        end
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,~,~] = tensorsvd(M,[1,2],[3,4],X);
    
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/2;
        C = C/max(C(:));
        
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/2;
        T = T/max(T(:));
        
        delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));  
        Cold = C;
        
        iter = iter + 1;
        if iter > maxiter
            disp(['Potts_CTM2 not converged at T = ' num2str(temp)]);
            break
        end            
    end       
end
 