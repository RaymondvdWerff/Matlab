function [Al,C,counter] = LeftOrthonormalize1(A,C,eta,temp)
    
    delta = 10;
    counter = 0;
    maxiter = 500;
    
    while delta > eta
        [~,C] = qrpos(C);
        C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

        Ac = ncon({C,A},{[-1,1],[1,-2,-3]});
        Cold = C;

        si = size(Ac);
        Ac = reshape(Ac,prod(si(1:2)),si(3));
        [Al,C] = qrpos(Ac);

        Al = reshape(Al,[si(1),si(2),si(3)]);
        C = reshape(C,si(3),si(3));
        lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
        C = C./lambda;
        delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));

        counter = counter+1;
        if counter > maxiter
            disp(['LeftOrthonormalize1 not converged at T = ' num2str(temp)]);
            break;
        end
    end
end