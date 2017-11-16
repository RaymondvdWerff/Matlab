function [Al,Ar,C1,C2,lambda,counter] = MPOFixedPoint1(T,A,m,eta)
    
    [Al,C1,~] = LeftOrthonormalize2(A,eye(size(A,2)),m,eta);
    [Ar,C2,~] = LeftOrthonormalize2(permute(A,[1,3,2]),eye(size(A,2)),m,eta);
    Ar = permute(Ar,[1,3,2]);
    C2 = permute(C2,[2,1]);
    
    Ml = ncon({Al,conj(Al),T},{[1,-4,-1],[2,-6,-3],[1,2,-5,-2]});
    siMl = size(Ml);
    Ml = reshape(Ml,prod(siMl(1:3)),prod(siMl(4:6)));
    Mr = ncon({Ar,conj(Ar),T},{[1,-1,-4],[2,-3,-6],[1,2,-2,-5]});
    siMr = size(Mr);
    Mr = reshape(Mr,prod(siMr(1:3)),prod(siMr(4:6)));
    
    [lambda,Fl] = LeadingEigenvector(Ml,rand(size(Ml,2),1),m,10^(-10));
    [~,Fr] = LeadingEigenvector(Mr,rand(size(Mr,2),1),m,10^(-10));
    Fl = reshape(Fl,[siMl(4),siMl(5),siMl(6)]);
    Fr = reshape(Fr,[siMr(4),siMr(5),siMr(6)]);
    
    C = C1*C2;
    Fl =  Fl./ncon({Fl,Fr,C,conj(C)},{[1,3,2],[4,3,5],[1,4],[2,5]});
    Ac = ncon({C1,A,C2},{[-2,1],[-1,1,2],[2,-3]});
    
    difference = ncon({Fl,Fr,Ac,T},{[2,3,-2],[4,5,-3],[1,2,4],[1,-1,3,5]})./lambda-ncon({Fl,Fr,C,Ar},{[1,2,4],[3,2,-3],[1,3],[-1,-2,4]});
    epsilon = sqrt(ncon({difference,conj(difference)},{[1,2,3],[1,2,3]})) 
    
    counter = 0;
    maxiter = 1000;
    
    while epsilon > eta
        M = ncon({Fl,Fr,T,C1,C2,inv(C1),inv(C2)},{[1,5,3],[2,6,4],[-4,-1,5,6],[1,-5],[-6,2],[-2,3],[4,-3]}); 
        siM = size(M);
        M = reshape(M,prod(siM(1:3)),prod(siM(4:6)));
        A = reshape(A,prod(siM(4:6)),1);
        [~,A] = LeadingEigenvector(M,A,m,epsilon/10);
        A = reshape(A,siM(4),siM(5),siM(6));
       
        [Al,C1,~] = LeftOrthonormalize2(A,C1,m,epsilon/10);
        [Ar,C2,~] = LeftOrthonormalize2(permute(A,[1,3,2]),permute(C2,[2,1]),m,epsilon/10);
        Ar = permute(Ar,[1,3,2]);
        C2 = permute(C2,[2,1]);

        Ml = ncon({Al,conj(Al),T},{[1,-4,-1],[2,-6,-3],[1,2,-5,-2]});
        siMl = size(Ml);
        Ml = reshape(Ml,prod(siMl(1:3)),prod(siMl(4:6)));
        Mr = ncon({Ar,conj(Ar),T},{[1,-1,-4],[2,-3,-6],[1,2,-2,-5]});
        siMr = size(Mr);
        Mr = reshape(Mr,prod(siMr(1:3)),prod(siMr(4:6)));

        [lambda,Fl] = LeadingEigenvector(Ml,reshape(Fl,size(Ml,2),1),m,10^(-10));
        [~,Fr] = LeadingEigenvector(Mr,reshape(Fl,size(Mr,2),1),m,10^(-10));
        Fl = reshape(Fl,[siMl(4),siMl(5),siMl(6)]);
        Fr = reshape(Fr,[siMr(4),siMr(5),siMr(6)]);

        C = C1*C2;
        Fl =  Fl./ncon({Fl,Fr,C,conj(C)},{[1,3,2],[4,3,5],[1,4],[2,5]});
        Ac = ncon({C1,A,C2},{[-2,1],[-1,1,2],[2,-3]});

        difference = ncon({Fl,Fr,Ac,T},{[2,3,-2],[4,5,-3],[1,2,4],[1,-1,3,5]})./lambda-ncon({Fl,Fr,C,Ar},{[1,2,4],[3,2,-3],[1,3],[-1,-2,4]});
        epsilon = sqrt(ncon({difference,conj(difference)},{[1,2,3],[1,2,3]})) 
        
        counter = counter + 1;
        if counter > maxiter
            disp('MPOFixedPoint1 not converged');
            break;
        end
    end
    
end