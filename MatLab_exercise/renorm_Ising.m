function [C,T,iter] = renorm_Ising(X,A,tol)
    
    rand(1);
    C = randn(X,X);
    C = (C + permute(C,[2,1]))/(2*max(abs(C(:))));
    T = randn(X,2,X);
    T = (T + permute(T,[3,2,1]))/(2*max(abs(T(:))));
    
    summ = 10;
    sv_old = randn(2*X);
    iter = 0;
    
    while summ > tol
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,sv,~] = tensorsvd(M,[1,2],[3,4],X);
        sv = sv/max(sv(:));
        
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/(2*max(abs(C(:))));
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/(2*max(abs(T(:))));
        
        summ = 0;
        for i = 1:X
            summ = summ + abs(sv(i,i)-sv_old(i,i));
        end
        
        sv_old = sv;  
        iter = iter + 1;
    end
end
    
    



