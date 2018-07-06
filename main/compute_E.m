function [E,iters,tictocs] = compute_E(Q,q,X,tol,maxiter,ts,func)
    
    delta_t = 1e-8;
    
    emptylist = zeros(numel(ts),1);    
    E = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end

    %Qsq = sqrtm(Q(q,ts(1),0));
    %A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    %[C,T] = beginmatrices(Qsq,A,X,1);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        tic
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        
        [C,T,iter1] = func(A,C,T,X,tol,maxiter,temp);
        E1 = log(compute_kappa(A,C,T));
        
        Qsq = sqrtm(Q(q,temp+delta_t,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        
        [C,T,iter2] = func(A,C,T,X,tol,maxiter,temp+delta_t);
        E2 = log(compute_kappa(A,C,T));
        
        E(t) = temp^2*(E2-E1)/delta_t;
        tictocs(t) = toc;
        iters(t) = iter1+iter2;      
    end
end
