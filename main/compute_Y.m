function [Y,iters,tictocs] = compute_Y(Q,q,X,tol,maxiter,ts,func)
    
    emptylist = zeros(numel(ts),1);    
    Y = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end

    for t = 1:numel(ts)
        
        temp = ts(t);
        disp(['temp = ' num2str(temp)]);
        
        Qsq = sqrtm(Q(q,temp,0));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C1,T1] = beginmatrices1(Qsq,A,X);
        [C2,T2] = beginmatrices2(Qsq,A,X);
        
        tic
        [C1,T1,iter1] = func(A,C1,T1,X,tol,maxiter,temp);
        [C2,T2,iter2] = func(A,C2,T2,X,tol,maxiter,temp);
        
        f1 = -temp*log(compute_kappa(A,C1,T1));
        f2 = -temp*log(compute_kappa(A,C2,T2));
        
        Y(t) = f2-f1;
        
        tictocs(t) = toc;
        iters(t) = iter1+iter2;      
    end
end

function k = compute_kappa2(A,C1,T1,C2,T2)
    
    env = ncon({C1,T1,C2,T2,C2,T1,C1,T1},{[1,2],[2,-1,3],[3,4],[4,-2,5],[5,6],[6,-3,7],[7,8],[8,-4,1]});
    k1 = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        
    k2 = ncon({C1,C1,C2,C2},{[1,2],[2,3],[3,4],[4,1]});

    k3 = ncon({C1,C2,T2,C2,C1,T1},{[1,2],[2,3],[3,4,5],[5,6],[6,7],[7,4,1]});

    k4 = ncon({C1,C1,T1,C2,C2,T1},{[1,2],[2,3],[3,4,5],[5,6],[6,7],[7,4,1]});
    
    k = k1*k2/(k3*k4);
    
end

function [C0,T0] = beginmatrices1(Qsq,A,X)
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

function [C0,T0] = beginmatrices2(Qsq,A,X)
    q = size(A,1);
    spin1_2D = zeros(q,q);spin1_2D(2,2) = 1;
    spin1_3D = zeros(q,q,q);spin1_3D(2,2,2)=1;
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