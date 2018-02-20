function [Xm,iters,tictocs] = compute_Xm(Q,q,X,tol,maxiter,ts,func)
    
    delta_h = 1e-8;
    
    [m1,iters1,tictocs1] = compute_m(Q,q,X,tol,maxiter,ts,func,0);
    [m2,iters2,tictocs2] = compute_m(Q,q,X,tol,maxiter,ts,func,delta_h);

    Xm = (m2-m1)/delta_h;
    iters = iters1+iters2;
    tictocs = tictocs1+tictocs2;
end

