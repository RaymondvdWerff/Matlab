function [Al,C] = LeftOrthonormalize6(A,tol,temp)

    X = size(A,1);
    q = size(A,2);

    %C2 = ncon({A,A},{[1,2,-2],[1,2,-1]});
    %C2 = arnoldi(A,A,reshape(C2,X^2,1),@mult1,min(X^2,15),tol);
    %C2 = reshape(C2,X,X);
    
    A2 = ncon({A,conj(A)},{[-3,1,-1],[-4,1,-2]});
    A2 = reshape(A2,X^2,X^2);

    [C2,~] = eigs(A2,1,'LM');
    C2 = reshape(C2,X,X);
    
    [U,s2,~] = svd(C2);
    C = U*sqrt(s2)*U';
    C = C./max(abs(C(:)));

    CA = ncon({C,A},{[-1,1],[1,-2,-3]});
    CA = reshape(CA,q*X,X);

    [U,s,V] = svd(CA,0);
    Al = U*V';
    Al = reshape(Al,[X,q,X]);
    C1 = V*s*V';
    C1 = C1./max(abs(C1(:)));
    
    delta = sqrt(ncon({C-C1,conj(C-C1)},{[1,2],[1,2]}));

    if delta > tol
            disp(['LeftOrthonormalize not converged at T = ' num2str(temp)]);
    end
end