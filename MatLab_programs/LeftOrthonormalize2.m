function [Al,C,counter] = LeftOrthonormalize2(A,C,m,eta,temp)

    counter = 0;

    [~,C] = qrpos(C);
    C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

    Ac = ncon({C,A},{[-2,1],[-1,1,-3]});
    Cold = C;

    si = size(Ac);
    Ac = reshape(Ac,prod(si(1:2)),si(3));
    [Al,C] = qrpos(Ac);    
    Al = reshape(Al,[si(1),si(2),si(3)]);

    lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
    C = C./lambda;
    delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));

    while delta > eta

        M = ncon({conj(Al),A},{[1,-3,-1],[1,-4,-2]});
        si = size(M);
        M = reshape(M,prod(si(1:2)),prod(si(3:4)));
        C = reshape(C,prod(si(3:4)),1);
        [~,C] = LeadingEigenvector(M,C,m,eta/10);
        
        C = reshape(C,si(3),si(4));

        [~,C] = qrpos(C);
        C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

        Ac = ncon({C,A},{[-2,1],[-1,1,-3]});
        Cold = C;

        si = size(Ac);
        Ac = reshape(Ac,prod(si(1:2)),si(3));
        [Al,C] = qrpos(Ac);
        Al = reshape(Al,[si(1),si(2),si(3)]);

        lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
        C = C./lambda;
        
        delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));
        counter = counter +1;
        if counter > maxiter
            disp(['LeftOrthonormalize2 not converged at T = ' num2str(temp)]);
            break;
        end
    end
end

