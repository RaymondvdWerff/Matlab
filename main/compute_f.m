function [f,iters,tictocs] = compute_f(Q,q,X,tol,maxiter,ts,func,h)
    
    emptylist = zeros(numel(ts),1);    
    f = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end

    Qsq = sqrtm(Q(q,ts(1),h));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    [C,T] = beginmatrices(Qsq,A,X);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        
        Qsq = sqrtm(Q(q,temp,h));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        f(t) = -temp*log(compute_kappa(A,C,T));

        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
