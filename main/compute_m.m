function [m,iters,tictocs] = compute_m(Q,q,X,tol,maxiter,ts,func,h)
    
    emptylist = zeros(1,numel(ts));    
    m = emptylist;
    %mx = emptylist;
    %my = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    %spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    
    Qsq = sqrtm(Q(q,ts(1),h));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    [C,T] = beginmatrices(Qsq,A,X);
    
    for t = 1:numel(ts)
        
        temp = ts(t)
        
        Qsq = sqrtm(Q(q,temp,h));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        %B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
        env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
        Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        mx = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]})/Z;
        my = ncon({env,By},{[1,2,3,4],[1,2,3,4]})/Z;
        
        %mx = collapse(C,T,Bx)/collapse(C,T,A);
        %my = collapse(C,T,By)/collapse(C,T,A);
        m(t) = sqrt(mx^2+my^2);
        
        %n = collapse(C,T,B)/collapse(C,T,A);
        %m(t) = (q*n-1)/(q-1);
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
