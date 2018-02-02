function [C,T,iter] = FPCM(A,C,T,X,tol,maxiter,temp)
    
    q = size(A,1);
    p = 12;
    for iter = 1:maxiter
        
        [Tl,C] = LeftOrthonormalize(T,tol,maxiter/10,temp,p);
        [~,s,~] = svd(C);s = s/max(s(:));
        
        opts.v0 = reshape(T,q*X^2,1);
        opts.tol = 1e-2;
        opts.p = p;
        [T,~] = eigs(@(x)mult2(Tl,A,x),q*X^2,1,'LM',opts);
        T = real(T); T = reshape(T,[X,q,X]);
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

function [Tl,C1] = LeftOrthonormalize(T,tol,maxiter,temp,p)

    X = size(T,1);
    q = size(T,2);
    
    for iter = 1:maxiter
        
        if iter == 1
            C2 = ncon({T,T},{[1,2,-2],[1,2,-1]});
            opts.v0 = reshape(C2,X^2,1);
            opts.tol = 1e-2;
            opts.p = p;
            [C2,~] = eigs(@(x)mult1(T,T,x),X^2,1,'LM',opts);
            C2 = reshape(C2,X,X);

            [U,s2,~] = svd(C2);
            C = U*sqrt(s2)*U';
            C = C./max(abs(C(:)));       
        else
            opts.v0 = reshape(C1,X^2,1);
            opts.tol = 1e-2;
            opts.p = p;
            [C,~] = eigs(@(x)mult1(T,Tl,x),X^2,1,'LM',opts);
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