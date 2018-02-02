function [c,iters,tictocs] = compute_cv(Q,q,X,tol,maxiter,ts,func,delta_b)
    
    emptylist = zeros(1,numel(ts));    
    c = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    Qsq = sqrtm(Q(q,ts(1)-delta_b,0));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    [C,T] = beginmatrices(Qsq,A,X);
    
    for t = 1:numel(ts)        
        tic
        
        b = 1/ts(t)-delta_b;
        Qsq = sqrtm(Q(q,1/b,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T,iter1] = func(A,C,T,X,tol,maxiter,1/b);
        k1 = compute_kappa(A,C,T);
        
        b = 1/ts(t);
        Qsq = sqrtm(Q(q,1/b,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T,iter2] = func(A,C,T,X,tol,maxiter,1/b);
        k2 = compute_kappa(A,C,T);
        
        b = 1/ts(t)+delta_b;
        Qsq = sqrtm(Q(q,1/b,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T,iter3] = func(A,C,T,X,tol,maxiter,1/b);
        k3 = compute_kappa(A,C,T);
        
        c(t) = ((log(k1)-2*log(k2)+log(k3))/delta_b^2)/ts(t)^2;
        iters(t) = iter1+iter2+iter3;
        tictocs(t) = toc;     
    end
end
