function [mx,my,mx2,my2,iters,tictocs] = compute_Xm3(Q,q,X,tol,maxiter,ts,func)
    
    emptylist = zeros(1,numel(ts));    
    Xm = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    spinx2_4D = zeros(q,q,q,q);for i=1:q; spinx2_4D(i,i,i,i)=sin(2*pi*(i-1)/q)^2; end
    spiny2_4D = zeros(q,q,q,q);for i=1:q; spiny2_4D(i,i,i,i)=cos(2*pi*(i-1)/q)^2; end
    
    %Qsq = sqrtm(Q(q,ts(1),h));
    %A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    %[C,T] = beginmatrices(Qsq,A,X);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        Bx2 = ncon({spinx2_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        By2 = ncon({spiny2_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        %C = rand(X);C = C+C'; T = rand(X,q,X); T = T + permute(T,[3,2,1]);
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
        env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
        Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        mx(t) = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]})/Z;
        my(t) = ncon({env,By},{[1,2,3,4],[1,2,3,4]})/Z;
        mx2(t) = ncon({env,Bx2},{[1,2,3,4],[1,2,3,4]})/Z;
        my2(t) = ncon({env,By2},{[1,2,3,4],[1,2,3,4]})/Z;
        
        %Xm(t) = mx2+my2-mx^2-my^2;
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
