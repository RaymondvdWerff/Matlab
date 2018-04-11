function [Xmx,Xmy,iters,tictocs] = compute_Xm(Q,q,X,tol,maxiter,ts,func)
    
    delta_h = 1e-8;
    
    [mx1,my1,iters1,tictocs1] = compute_m(Q,q,X,tol,maxiter,ts,func,0);
    [mx2,my2,iters2,tictocs2] = compute_m(Q,q,X,tol,maxiter,ts,func,delta_h);

    Xmx = (mx2-mx1)/delta_h;
    Xmy = (my2-my1)/delta_h;
    iters = iters1+iters2;
    tictocs = tictocs1+tictocs2;
end

