function [Al,C] = LeftOrthonormalize3(A,C,eta,temp)
    
    iter = 0;
    maxiter = 500;
    D = size(A,1);
    
    [~,C] = qrpos(C);
    C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

    Ac = ncon({C,A},{[-1,1],[1,-2,-3]});
    Cold = C;

    si = size(Ac);
    Ac = reshape(Ac,prod(si(1:2)),si(3));
    [Al,C] = qrpos(Ac);    
    Al = reshape(Al,[si(1),si(2),si(3)]);
    
    lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
    C = C./lambda;
    
    delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));

    while delta > eta
        
        C = reshape(C,D^2,1);
        opts.v0 = C;
        opts.isreal = false;
        
        [C,~] = eigs(@(x)leadingvec(x,A,Al),D^2,1,'LM',opts);
        C = reshape(C,D,D);
        
        [~,C] = qrpos(C);
        C = C./sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));

        Ac = ncon({C,A},{[-1,1],[1,-2,-3]});
        Cold = C;

        si = size(Ac);
        Ac = reshape(Ac,prod(si(1:2)),si(3));
        [Al,C] = qrpos(Ac);
        
        Al = reshape(Al,[si(1),si(2),si(3)]);

        lambda = sqrt(ncon({C,conj(C)},{[1,2],[1,2]}));
        C = C./lambda;
        
        delta = sqrt(ncon({C-Cold,conj(C-Cold)},{[1,2],[1,2]}));
        
        iter = iter +1;
        if iter > maxiter
            disp(['LeftOrthonormalize3 not converged at T = ' num2str(temp)]);
            break;
        end
    end
end

function y = leadingvec(x,A,Al)
    D = size(A,1);
    x = reshape(x,D,D);
    y = ncon({conj(Al),x,A},{[1,3,-1],[1,2],[2,3,-2]});
    y = reshape(y,D^2,1);
end
