function [c,tictocs] = compute_cv(Q,q,X,tol,maxiter,ts,func,delta_b)
    
    emptylist = zeros(1,numel(ts));    
    c = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    Qsq = sqrtm(Q(q,ts(1)-delta_b,0));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    [C,T] = beginmatrices(Qsq,A,X);
    
    for t = 1:numel(ts)        
        tic
        
        b = 1/ts(t)-delta_b;
        Qsq = sqrtm(Q(q,1/b,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T,~] = func(A,C,T,X,tol,maxiter,1/b);
        k1 = compute_kappa(A,C,T);
        
        b = 1/ts(t);
        Qsq = sqrtm(Q(q,1/b,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T,~] = func(A,C,T,X,tol,maxiter,1/b);
        k2 = compute_kappa(A,C,T);
        
        b = 1/ts(t)+delta_b;
        Qsq = sqrtm(Q(q,1/b,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T,~] = func(A,C,T,X,tol,maxiter,1/b);
        k3 = compute_kappa(A,C,T);
        
        c(t) = ((log(k1)-2*log(k2)+log(k3))/delta_b^2)/ts(t)^2;

        tictocs(t) = toc;     
    end
end
    
function [C0,T0] = beginmatrices(Qsq,A,X)
    q = size(A,1);
    spin1_2D = zeros(q,q);spin1_2D(1,1) = 1;
    spin1_3D = zeros(q,q,q);spin1_3D(1,1,1)=1;
    C0 = ncon({spin1_2D,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T0 = ncon({spin1_3D,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    
    while size(T0,1) < X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TA = ncon({T0,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TA},{[-1,1,2],[1,2,-2,-3,-4]});
        C0 = reshape(M,q*size(T0,1),q*size(T0,1));
        T0 = reshape(TA,[q*size(T0,1),q,q*size(T0,1)]);
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
    
    if size(T0,1) > X
        [U,~,~] = svd(C0);
        U_til = U(:,1:X);
        C0 = ncon({C0,U_til,U_til},{[1,2],[1,-1],[2,-2]});
        T0 = ncon({T0,U_til,U_til},{[1,-2,2],[1,-1],[2,-3]});
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
end