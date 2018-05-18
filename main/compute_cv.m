function [cv,iters,tictocs] = compute_cv2(Q,q,X,tol,maxiter,ts,func)
    
    delta_t = 1e-5;
    
    emptylist = zeros(1,numel(ts));    
    cv = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    
    for t = 1:numel(ts)        
        tic
        
        disp(['temp = ' num2str(ts(t))]);
        
        temp = ts(t)-delta_t;
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        [C,T,iter1] = func(A,C,T,X,tol,maxiter,temp);
        f1 = -temp*log(compute_kappa(A,C,T));
        
        temp = ts(t);
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        [C,T,iter2] = func(A,C,T,X,tol,maxiter,temp);
        f2 = -temp*log(compute_kappa(A,C,T));
        
        temp = ts(t)+delta_t;
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        [C,T,iter3] = func(A,C,T,X,tol,maxiter,temp);
        f3 = -temp*log(compute_kappa(A,C,T));
        
        cv(t) = -ts(t)*(f1-2*f2+f3)/delta_t^2;
        iters(t) = iter1+iter2+iter3;
        tictocs(t) = toc;     
    end
end
