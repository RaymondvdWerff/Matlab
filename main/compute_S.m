function [S,iters,tictocs] = compute_S(Q,q,X,tol,maxiter,ts,func)
    
    emptylist = zeros(1,numel(ts));    
    S = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    
    Qsq = sqrtm(Q(q,ts(1),0));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    [C,T] = beginmatrices(Qsq,A,X,1);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        %C = rand(X);C = C+C'; T = rand(X,q,X); T = T + permute(T,[3,2,1]);
        %[C,T] = beginmatrices(Qsq,A,X,0);
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        v = eig(C);v = v/max(v);
        
        S(t) = sum(-v.^4.*log(v.^4));
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
