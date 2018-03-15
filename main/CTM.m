function [C,T,iter] = CTM(A,C,T,X,tol,maxiter,temp)

    for iter = 1:maxiter
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,s,~] = tensorsvd(M,[1,2],[3,4],X);
        s = s/max(s(:));
    
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/max(C(:));
        
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/max(T(:));
        
        if iter > 1
            delta = sum(sum(abs(s-sold)));
            if delta < tol
                break;
            end
        end  
        sold = s;           
    end
    
    if (iter == maxiter) && (maxiter > 10)
        disp(['CTM not converged at T = ' num2str(temp)]);
    end 
end
 