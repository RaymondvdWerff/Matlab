function [C,T,iter] = FPCM2(A,C,T,X,tol,maxiter,temp)
    
    q = size(A,1);
    
    for iter = 1:maxiter
        
        [Tl,C] = LeftOrthonormalize(T,tol,maxiter/10,temp);
        [~,s,~] = svd(C);s = s/max(s(:));

        T = leadingvec(Tl,A,T,@mult2);
        T = reshape(T,[X,q,X]);
        T = T + permute(T,[3,2,1]);

        if iter > 1
            delta = sum(sum(abs(s-sold)));
            if delta < tol
                break;
            end
        end  
        sold = s;       
    end
    
    if iter == maxiter
        disp(['FPCM not converged at T = ' num2str(temp)]);
    end 
end

function [Tl,C1] = LeftOrthonormalize(T,tol,maxiter,temp)

    X = size(T,1);
    q = size(T,2);
    
    for iter = 1:maxiter
        
        if iter == 1
            C2 = ncon({T,T},{[1,2,-2],[1,2,-1]});
            C2 = leadingvec(T,T,C2,@mult1);
            C2 = reshape(C2,X,X);

            [U,s2,~] = svd(C2);
            C = U*sqrt(s2)*U';
            C = C./max(abs(C(:)));       
        else
            C = leadingvec(T,Tl,C1,@mult1);
            C = reshape(C,X,X);

            [~,s,V] = svd(C);
            C = V*s*V';
            C = C/max(abs(C(:)));
        end

        CT = ncon({C,T},{[-1,1],[1,-2,-3]});
        CT = reshape(CT,q*X,X);

        [U,s,V] = svd(CT,0);
        Tl = U*V';
        Tl = reshape(Tl,[X,q,X]);
        C1 = V*s*V';
        C1 = C1/max(abs(C1(:)));
        delta = sqrt(ncon({abs(C-C1),abs(C-C1)},{[1,2],[1,2]}));
        
        if delta < tol
            break
        end
    end
    if iter == maxiter
        disp(['LeftOrthonormalize not converged at T = ' num2str(temp)]);
    end
end

function y = leadingvec(A,B,C,func)
    
    n = numel(C);
    m = 15;
    V = zeros(n,m);
    T = zeros(m,m);

    v = reshape(C,n,1);
    v = v/norm(v);
    V(:,1) = v;
    w_prime = func(A,B,V(:,1));
    T(1,1) = w_prime'*V(:,1);
    w = w_prime - T(1,1)*V(:,1);

    for j = 2:m
        T(j-1,j) = norm(w);
        if T(j-1,j) == 0
            disp('special case');
        end
        T(j,j-1) = T(j-1,j);
        V(:,j) = w/T(j-1,j);
        w_prime = func(A,B,V(:,j));
        T(j,j) = w_prime'*V(:,j);
        w = w_prime - T(j,j)*V(:,j) - T(j-1,j)*V(:,j-1);
        for i = 1:j-1
           w = w - w'*V(:,i)*V(:,i); 
        end
    end

    Qtot = eye(m);
    for i = 1:5
        [Q,R] = qr(T);
        Qtot = Qtot*Q;
        T = R*Q;
    end
    
    y = V*Qtot(:,1);
end

function y = mult1(M1,M2,C)
    si = size(M1);
    C = reshape(C,si(1),si(1));
    y = ncon({M1,M2,C},{[1,3,-2],[2,3,-1],[2,1]});
    y = reshape(y,si(1)^2,1);
end

function y = mult2(M1,M2,T)
    si = size(M1);
    T = reshape(T,[si(1),si(2),si(3)]);
    y = ncon({T,M1,M2,M1},{[4,1,2],[2,3,-3],[3,5,1,-2],[4,5,-1]});
    y = reshape(y,prod(si),1);
end