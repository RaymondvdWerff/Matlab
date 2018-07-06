function [f,iters,tictocs] = compute_f(Q,q,X,tol,maxiter,ts,func)
    emptylist = zeros(numel(ts),1);    
    f = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        f(t) = -temp*log(compute_kappa(A,C,T));

        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
