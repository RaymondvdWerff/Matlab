function [mx,my,corr,iters,tictocs] = compute_m_and_corr(Q,q,X,tol,maxiter,ts,func,h)
    
    emptylist = zeros(1,numel(ts));    
    mx = emptylist;
    my = emptylist;
    corr = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
   
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,h));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
        env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
        Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        
        mx(t) = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]})/Z;
        my(t) = ncon({env,By},{[1,2,3,4],[1,2,3,4]})/Z;
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