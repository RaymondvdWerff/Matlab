function [iters,ticks,m] = converge_m_CTM2(C,T,A0,A1,q,X,maxiter,tol)

    opts.isreal = false;
    tic
    
    for iter = 1:maxiter
        
        if iter == 1
            A = A0;
        end
        if iter > 1
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