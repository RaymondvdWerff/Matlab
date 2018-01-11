function [corr,iters,tictocs] = compute_corr2(Q,q,X,tol,maxiter,ts,func)
    
    emptylist = zeros(1,numel(ts));    
    corr = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    
    C = rand(X); C = C + C';
    T = rand(X,q,X); T = T + permute(T,[3,2,1]);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        Qsq = sqrtm(Q(q,temp));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        corr(t) = corrlen(T);
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end

function y = corrlen(T)
    X = size(T,1);
    [~,vals] = eigs(@(x)mult1(T,T,x),X^2,2,'LM');
    y = 1/abs(log(abs(vals(1,1))/abs(vals(2,2))));
end

function y = mult1(M1,M2,C)
    si = size(M1);
    C = reshape(C,si(1),si(1));
    y = ncon({M1,M2,C},{[1,3,-2],[2,3,-1],[2,1]});
    y = reshape(y,si(1)^2,1);
end
