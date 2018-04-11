function [mx,my,iters,tictocs] = compute_m(Q,q,X,tol,maxiter,ts,func,h)
    
    emptylist = zeros(1,numel(ts));    
    mx = emptylist;
    my = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    %spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    %spin1 = zeros(q,q,q,q);spin1(1,1,1,1) = 1;spin1(2,2,2,2) = -1;
    %spin2 = zeros(q,q,q,q);spin2(3,3,3,3) = 1;spin2(4,4,4,4) = -1;
    
    %Qsq = sqrtm(Q(q,ts(1),0));
    %A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    %[C,T] = beginmatrices(Qsq,A,X,1);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,h));Qsq = real(Qsq);
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        %B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        %Bx = ncon({spin1,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        %Bz = ncon({spin2,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
        tic
        [C,T] = beginmatrices(Qsq,A,X,1);
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
        env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
        Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        
        mx(t) = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]})/Z;
        my(t) = ncon({env,By},{[1,2,3,4],[1,2,3,4]})/Z;
        %n = ncon({env,B},{[1,2,3,4],[1,2,3,4]})/Z;
        %m(t) = (q*n-1)/(q-1);
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
