function [C1,C2,T1,T2,i] = renorm_Ising_triang(Qsq,X,A,tol)

    maxsteps=20000;
    minsteps=4;
    
    q = size(A,2);
    [C1,C2,T1,T2] = beginmatrices(Qsq,q);

    sold = zeros(q,1);

    for i=1:maxsteps

        N1 = ncon({T1,A},{[-1,-5,-2,1],[1,-3,-4,-6,-7,-8]});
        M1 = ncon({C1,T2,N1},{[1,2],[-1,1,3,4],[2,3,4,-2,-4,-5,-6,-3]});
     
        N2 = ncon({T2,A},{[-1,-4,-5,1],[1,-2,-3,-6,-7,-8]});
        M2 = ncon({C2,T1,N2},{[1,2,3],[-1,1,-2,4],[2,3,4,-4,-5,-6,-7,-3]});
        
        si1 = size(M1);
        M_SVD1 = reshape(M1,prod(si1(1:3)),prod(si1(4:6)));
        M_SVD1 = (M_SVD1+M_SVD1')/2;
        D = min(X,size(M_SVD1,1));
        
        [U1,s1,~] = svd(M_SVD1,'econ');
        
        U1 = U1(:,1:D);
        U1 = reshape(U1,[si1(1),si1(2),si1(3),D]);
        s1 = diag(s1);
        s1 = s1(1:D);
        snew = s1/max(s1);
    
        C1 = ncon({M1,conj(U1),U1},{[1,2,3,4,5,6],[1,2,3,-1],[4,5,6,-2]});
        C1 = C1+permute(C1,[2 1]);
        C1 = C1./max(abs(C1(:)));
        
        T1 = ncon({conj(U1),N1,U1},{[1,2,3,-1],[1,2,3,-3,4,5,6,-4],[4,5,6,-2]}); 
        T1 = T1+permute(T1,[2 1 3 4]);
        T1 = T1./max(abs(T1(:)));
        
        C2 = ncon({M2,conj(U1),U1},{[1,2,3,4,5,6,-3],[1,2,3,-1],[4,5,6,-2]});
        C2 = C2+permute(C2,[2 1 3]);
        C2 = C2./max(abs(C2(:)));
        
        T2 = ncon({conj(U1),N2,U1},{[1,2,3,-1],[1,2,3,4,5,6,-3,-4],[4,5,6,-2]}); 
        T2 = T2+permute(T2,[2 1 3 4]);
        T2 = T2./max(abs(T2(:)));
        
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

function [C1_0,C2_0,T1_0,T2_0] = beginmatrices(Qsq,q)
    C1_0 = zeros(q,q);
    C1_0(1,1) = 1;
    C1_0 = ncon({C1_0,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    
    C2_0 = zeros(q,q,q);
    C2_0(1,1,1) = 1;
    C2_0 = ncon({C2_0,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    
    T1_0 = zeros(q,q,q,q);
    T1_0(1,1,1,1) = 1;
    T1_0 = ncon({T1_0,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});

    T2_0 = T1_0;
end





