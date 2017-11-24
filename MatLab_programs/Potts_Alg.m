function [T1,C,iter] = Potts_Alg(T1,C_prime,A,A1,q,X,tol,temp)
    
    delta = tol + 1;
    iter = 0;
    maxiter = 700;
    opts.isreal = false;
       
    while delta > tol
        
        if iter > 0
            A = A1;
        end
        
        [T1l,C] = LeftOrthonormalize3(T1,C_prime,tol,temp);
        T1l = real(T1l);C = real(C);
        
        [T2,~] = eigs(@(x)leadingvec(x,T1l,A),q*X^2,1,'LM');     
        T2 = reshape(T2,[X,q,X]);
        T2 = T2 + permute(T2,[3,2,1]);
        
        [T2r,C_prime] = LeftOrthonormalize3(permute(T2,[3,2,1]),C,tol,temp);
        T2r = real(T2r);C_prime = real(C_prime);
        T2r = permute(T2r,[3,2,1]);
        C_prime = permute(C_prime,[2,1]);
        
        [T1,~] = eigs(@(x)leadingvec(x,permute(T2r,[3,2,1]),A),q*X^2,1,'LM');
        T1 = reshape(T1,[X,q,X]);
        T1 = T1 + permute(T1,[3,2,1]);
        
        delta = sqrt(ncon({C-C_prime,conj(C-C_prime)},{[1,2],[1,2]})); 
        
        iter = iter + 1;
        if iter > maxiter
            disp(['Potts_Alg not converged at T = ' num2str(temp)]);
            break
        end
    end
end

function y = leadingvec(x,M1,M2)
    si = size(M1);
    x = reshape(x,[si(1),si(2),si(3)]);
    y = ncon({x,M1,M2,M1},{[4,1,2],[2,3,-3],[3,5,1,-2],[4,5,-1]});
    y = reshape(y,prod(si),1);
end
