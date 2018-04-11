function [mx,mz,iters] = compute_mth(Q,q,X,tol,maxiter,ts,hs,func)
    
    emptylist = zeros(numel(ts),numel(hs));    
    mx = emptylist;
    mz = emptylist;
    iters = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    spin1 = zeros(q,q,q,q);spin1(1,1,1,1) = 1;spin1(2,2,2,2) = -1;
    spin2 = zeros(q,q,q,q);spin2(3,3,3,3) = 1;spin2(4,4,4,4) = -1;
    
    for t = 1:numel(ts)
        for h = 1:numel(hs)

            temp = ts(t);
            magf = hs(h);
            
            disp(['temp = ' num2str(temp) ' ,magf = ' num2str(magf)]);

            Qsq = sqrtm(Q(q,magf,temp));Qsq = real(Qsq);
            A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
            Bx = ncon({spin1,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
            Bz = ncon({spin2,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});

            [C,T] = beginmatrices(Qsq,A,X,1);

            tic
            [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);

            env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
            env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
            Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
            mx(t,h) = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]})/Z;
            mz(t,h) = ncon({env,Bz},{[1,2,3,4],[1,2,3,4]})/Z;
            iters(t,h) = iter;
        end
    end
end
