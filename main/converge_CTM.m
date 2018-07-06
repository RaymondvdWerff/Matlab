function [f,m,sv,iter,tictocs] = converge_CTM(Q,q,X,maxiter,temp)
    
    f = zeros(1,maxiter);
    m = zeros(1,maxiter);
    sv = zeros(1,maxiter);
    tictocs = zeros(1,maxiter);
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    
    Qsq = sqrtm(Q(q,temp,0));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    [C,T] = beginmatrices(Qsq,A,X,1);
       
    for iter = 1:maxiter
        
        tic
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,s,~] = tensorsvd(M,[1,2],[3,4],X);s = s/max(s(:));
    
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/max(C(:));
        
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/max(T(:));
        
        if iter == 1
            tictocs(iter) = toc;
        else
            tictocs(iter) = tictocs(iter-1) + toc;
        end
        
        [k,M] = compute_kappa_m(A,Bx,By,C,T);
        f(iter) = -temp*log(k);
        m(iter) = M;
        
        if iter > 1
            sv(iter) = sum(sum(abs(s-sold)));
        end
        if iter > 10
            if (abs(m(iter)-m(iter-10)) < 1e-15)
                break
            end
        end
        sold = s;
    end
    if iter == maxiter
        disp('converge_CTM not converged');
    end
end
