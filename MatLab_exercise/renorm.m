function [C,T,iter] = renorm(q,b,X,A,tol)    
        
    Qsq = sqrtm(Q(q,b));
    [C0,T0] = beginmatrices(q,1);
    
    C = ncon({C0,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T = ncon({T0,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
  
    sv = randn(q*X);
    summ = 10;    
    iter = 0;
    dim = q;
    
    %{
    if X > q
        while dim < X
            dim = dim*q;
            C_new = (C + permute(C,[2,1]))/2;
            C = C_new/max(C_new(:));
            T_new = (T + permute(T,[3,2,1]))/2;
            T = T_new/max(T_new(:));
            C_old = ncon({C,T,T,A},{[1,2],[2,3,-3],[1,4,-1],[4,3,-2,-4]});
            C = reshape(C_old,[dim,dim]);
            T_old = ncon({T,A},{[-1,1,-4],[-2,-3,-5,1]});
            T = reshape(T_old,[dim,q,dim]);
        end
        
        [U,~,~] = svd(C);
        
        U_til = U(:,1:X);
        U_til = reshape(U_til,[dim/q,q,X]);
        
        C_new = ncon({U_til,C_old,U_til},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C_new = (C_new + permute(C_new,[2,1]))/2;
        C = C_new/max(C_new(:));
        
        T_new = ncon({U_til,T_old,U_til},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T_new = (T_new + permute(T_new,[3,2,1]))/2;
        T = T_new/max(T_new(:));
    end
    
    if X < q
        M = ncon({C,T,T,A},{[1,2],[2,3,-3],[1,4,-1],[4,3,-2,-4]});
        M_SVD = reshape(M,[q^2,q^2]);
        N = ncon({T,A},{[-1,1,-4],[-2,-3,-5,1]});
        
        [U,~,~] = svd(M_SVD);
        
        U_til = U(:,1:X);
        U_til = reshape(U_til,[q,q,X]);
    
        C_new = ncon({U_til,M,U_til},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C_new = (C_new + permute(C_new,[2,1]))/2;
        C = C_new/max(C_new(:));
        
        T_new = ncon({U_til,N,U_til},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T_new = (T_new + permute(T_new,[3,2,1]))/2;
        T = T_new/max(T_new(:));
    end
    %}
    
    C = rand(X,X);
    T = 
    
    while summ > tol
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,sv_new,~] = tensorsvd(M,[1,2],[3,4],X);
        sv_new = sv_new/max(sv_new(:));
    
        C_new = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C_new = (C_new + permute(C_new,[2,1]))/2;
        C = C_new/max(C_new(:));
        
        T_new = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T_new = (T_new + permute(T_new,[3,2,1]))/2;
        T = T_new/max(T_new(:));
        
        summ = 0;
        for k = 1:sqrt(numel(sv_new))
            summ = summ + abs(sv(k,k)-sv_new(k,k));
        end
        
        sv = sv_new;
        iter = iter + 1;
            
    end       
    
end


function [C,T] = beginmatrices(q,s)
    C = zeros(q,q);
    C(s,s) = 1;
    
    T = zeros(q,q,q);
    T(s,s,s) = 1;        
end
       