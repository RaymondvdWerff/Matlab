function [Al,lambda,counter] = MPOFixedPoint3(T,A,m,eta)
    
    delta = 10;
    C1_prime = eye(size(A,2)); 
    counter = 0;
    maxiter  = 100;
    
    while delta > eta
       
        [Al,C1,~] = LeftOrthonormalize2(A,C1_prime,m,eta);
        C1 = C1./sqrt(ncon({C1,conj(C1)},{[1,2],[1,2]}));
        
        M1 = ncon({Al,conj(Al),T},{[1,-4,-1],[2,-6,-3],[1,2,-5,-2]});
        siM1 = size(M1);
        M1 = reshape(M1,prod(siM1(1:3)),prod(siM1(4:6)));
        [lambda,B] = LeadingEigenvector(M1,rand(size(M1,2),1),m,eta);
        B = reshape(B,[siM1(4),siM1(5),siM1(6)]);
        
        [Br,C1_prime] = LeftOrthonormalize2(permute(B,[1,3,2]),C1,m,eta);
        Br = permute(Br,[1,3,2]);
        C1_prime = permute(C1_prime,[2,1]);
        C1_prime = C1_prime./sqrt(ncon({C1_prime,conj(C1_prime)},{[1,2],[1,2]}));
        
        M2 = ncon({conj(Br),Br,T},{[-6,2,-3],[-5,1,-2],[-4,-1,1,2]});
        siM2 = size(M2);
        M2 = reshape(M2,prod(siM2(1:3)),prod(siM2(4:6)));
        [lambda_prime,A] = LeadingEigenvector(M2,rand(size(M2,2),1),m,eta);
        A = reshape(A,[siM2(4),siM2(5),siM2(6)]);
        
        delta = sqrt(ncon({C1-C1_prime,conj(C1-C1_prime)},{[1,2],[1,2]}))
        
        counter = counter + 1;
        if counter > maxiter
            disp('MPOFixedPoint3 not converged');
            break;
        end     
    end
end