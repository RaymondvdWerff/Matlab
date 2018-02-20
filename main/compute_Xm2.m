function [Xm,Xmx,Xmy,iters,tictocs] = compute_Xm2(Q,q,X,tol,maxiter,ts,func)
    
    delta_h = 1e-8;
    
    [mx1,my1,iters1,tictocs1] = compute_m2(Q,q,X,tol,maxiter,ts,func,0);
    [mx2,my2,iters2,tictocs2] = compute_m2(Q,q,X,tol,maxiter,ts,func,delta_h);

    Xm = (sqrt(mx2.^2+my2.^2)-sqrt(mx1.^2+my1.^2))/delta_h;
    Xmx = (mx2-mx1)/delta_h;
    Xmy = (my2-my1)/delta_h;
    iters = iters1+iters2;
    tictocs = tictocs1+tictocs2;
end

         