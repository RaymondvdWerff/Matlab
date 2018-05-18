function [evs,iters,tictocs] = compute_evs(Q,q,X,tol,maxiter,ts,func,numevs)
    
    evs = zeros(numel(ts),numevs);
    iters = zeros(numel(ts),1);
    tictocs = zeros(numel(ts),1);
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X,1);
        
        tic
        [~,T,iter] = func(A,C,T,X,tol,maxiter,temp);
        
        [~,ev] = eigs(@(x)mult1(T,T,x),X^2,numevs,'LM');
        ev = sort(diag(ev),'descend');ev = ev/ev(1);
        evs(t,:) = ev;
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end

function y = mult1(M1,M2,C)
    si = size(M1);
    C = reshape(C,si(1),si(1));
    y = ncon({M1,M2,C},{[1,3,-2],[2,3,-1],[2,1]});
    y = reshape(y,si(1)^2,1);
end
