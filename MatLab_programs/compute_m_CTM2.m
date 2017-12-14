function [ms,iters,tictocs] = compute_m_CTM2(q,X,tol,maxiter,ts)
    
    emptylist = zeros(1,numel(ts));    
    ms = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    spin1 = zeros(q,q,q,q);spin1(1,1,1,1)=1;
    
    for i = 1:numel(ts)
        
        t = ts(i);
        Qsq = sqrtm(ones(q)+(exp(1/t)-1)*eye(q));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        B = ncon({spin1,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X);
        
        tic;
        for iter = 1:maxiter

            N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
            M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});

            [U,s,~] = tensorsvd(M,[1,2],[3,4],X);
            s = s./max(s(:));
            
            C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
            C = (C + permute(C,[2,1]))/2;
            C = C/max(C(:));

            T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
            T = (T + permute(T,[3,2,1]))/2;
            T = T/max(T(:));

            n = collapse(C,T,B)/collapse(C,T,A);
            m = (q*n-1)/(q-1);

            if iter > 1
                delta = sum(sum(abs(s-sold)));
                if delta < tol
                    break;
                end
            end
            mold = m;
            sold = s;
        end
        
        if iter == maxiter
            disp(['compute_m_CMT2 not converged at T = ' num2str(t)]);
        end

        tictocs(i) = toc;
        ms(i) = m;
        iters(i) = iter;        
    end
end
    
function [C0,T0] = beginmatrices(Qsq,A,X)
    q = size(A,1);
    delta_2D = zeros(q,q);delta_2D(1,1)=1;
    delta_3D = zeros(q,q,q);delta_3D(1,1,1)=1;
    C0 = ncon({delta_2D,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T0 = ncon({delta_3D,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    
    while size(T0,1) < X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TA = ncon({T0,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TA},{[-1,1,2],[1,2,-2,-3,-4]});
        C0 = reshape(M,q*size(T0,1),q*size(T0,1));
        T0 = reshape(TA,[q*size(T0,1),q,q*size(T0,1)]);
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
    
    if size(T0,1) > X
        [U,~,~] = svd(C0);
        U_til = U(:,1:X);
        C0 = ncon({C0,U_til,U_til},{[1,2],[1,-1],[2,-2]});
        T0 = ncon({T0,U_til,U_til},{[1,-2,2],[1,-1],[2,-3]});
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
end

function v = arnoldi(M1,M2,b,m,tol)
    
    H = zeros(m+1,m);
    V = zeros(size(b,1),m+1);
    V(:,1) = b/norm(b);
    
    for j = 1:m
        w = mult(M1,M2,V(:,j));
        for i = 1:j
            H(i,j) = w'*V(:,i);
            w = w - H(i,j)*V(:,i);
        end
        
        H(j+1,j) = norm(w);
        
        if H(j+1,j) < tol
            H = H(1:j,1:j);
            break
        end
        
        V(:,j+1) = w./H(j+1,j);
    end
    
    H = H(1:j,1:j);
    V = V(:,1:j);
    [vector,~] = eigs(H,1,'LM');
    v = V*vector;
    
end

function y = mult(M1,M2,x)
    si = size(M1);
    x = reshape(x,[si(1),si(2),si(3)]);
    y = ncon({x,M1,M2,M1},{[2,1,4],[-3,3,2],[3,5,-2,1],[-1,5,4]});
    y = reshape(y,prod(si),1);
end