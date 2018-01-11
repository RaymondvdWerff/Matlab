function [C,T,i] = renorm_Ising2(Qsq,X,A,tol)

    maxsteps=20000;
    minsteps=4;
   
    q = size(A,2);
    [C,T] = beginmatrices(Qsq,q);

    sold = zeros(q,1);

    for i=1:maxsteps
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[-1,3,1],[2,3,-2,-3,-4]});
     
        si = size(M);
        M_SVD = reshape(M,prod(si(1:2)),prod(si(3:4)));
        M_SVD = (M_SVD+M_SVD')/2;
        D = min(X,size(M_SVD,1));

        [U,s,V] = svd(M_SVD,'econ');
        
        U = U (:,1:D);
        U = reshape(U,[si(1) si(2) D]);
        s = diag(s);
        s = s(1:D);
        snew = s/max(s);
    
        C = ncon({M,conj(U),U},{[1,2,3,4],[1,2,-1],[3,4,-2]});
        C = C+permute(C,[2 1]);
        C = C./max(abs(C(:)));
        
        T = ncon({U,N,conj(U)},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]}); 
        T=T+permute(T,[3 2 1]);
        T = T./max(abs(T(:)));
        
        if numel(sold) == numel(snew)
            diffs =  norm(snew-sold);
        else
            diffs = inf;
        end

        if (diffs<tol)  && (i>minsteps) 
            break;
        end
        sold = snew;
    end

end 

function [C0,T0] = beginmatrices(Qsq,q)
    C0 = zeros(q,q);
    C0(1,1) = 1;
    C0 = ncon({C0,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    
    T0 = zeros(q,q,q);
    T0(1,1) = 1;
    T0 = ncon({T0,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
end





