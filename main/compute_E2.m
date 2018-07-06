kfunction [E,iters,tictocs] = compute_E2(Q,q,X,tol,maxiter,ts,func,h)
    
    emptylist = zeros(1,numel(ts));    
    E = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
        
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,h));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        B = compute_E_4D(Qsq,q);
        
        tic
        [C,T] = beginmatrices(Qsq,A,X,1);
        
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        tictocs(t) = toc;
        iters(t) = iter;
        
        env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
        env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
        Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});

        E(t) = ncon({env,B},{[1,2,3,4],[1,2,3,4]})/Z;
    end
end

function [energy] = compute_E_4D(Qsq,q)
    energy = zeros(q,q,q,q);
    for i = 1:q
        for j = 1:q
            for k = 1:q
                for l = 1:q
                    for n = 1:q
                        energy(i,j,k,l) = energy(i,j,k,l)+(-cos(2*pi*(i-n)/q)-cos(2*pi*(j-n)/q)-cos(2*pi*(k-n)/q)-cos(2*pi*(l-n)/q))*Qsq(i,n)*Qsq(j,n)*Qsq(k,n)*Qsq(l,n)/2;
                    end
                end
            end
        end
    end
end
