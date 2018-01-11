function [m,iters,tictocs,imarkers,tmarkers] = converge_m_CTM(Q,q,X,tol,maxiter,temp,tols)
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    %spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    
    Qsq = sqrtm(Q(q,temp));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    %B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    
    [C,T] = beginmatrices(Qsq,A,X);
    i = 1;
    tictocs = [0];
       
    for iter = 1:maxiter
        tic;
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,s,~] = tensorsvd(M,[1,2],[3,4],X);s = s/max(s(:));
    
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/max(C(:));
        
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/max(T(:));
        
        tictocs(iter) = tictocs(end) + toc;
        
        mx = collapse(C,T,Bx)/collapse(C,T,A);
        my = collapse(C,T,By)/collapse(C,T,A);
        m(iter) = sqrt(mx^2+my^2);
        %n = collapse(C,T,B)/collapse(C,T,A);
        %m(iter) = (q*n-1)/(q-1);
        
        
        if iter > 1
            %m(iter) = sum(sum(abs(s-sold)));
            if i < numel(tols)+1
                if sum(sum(abs(s-sold))) < tols(i)
                    imarkers(i) = iter;
                    tmarkers(i) = tictocs(iter);
                    i = i+1;
                end
            end
            if abs(m(iter)-m(iter-1)) < tol
                break;
            end
        end
        sold = s;
    end
    
    if iter == maxiter
        disp('converge_m_CTM not converged');
    end 
    iters = 1:iter;
end

function [C0,T0] = beginmatrices(Qsq,A,X)
    q = size(A,1);
    spin1_2D = zeros(q,q);spin1_2D(1,1)=1;
    spin1_3D = zeros(q,q,q);spin1_3D(1,1,1)=1;
    spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    C0 = ncon({spin1_2D,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T0 = ncon({spin1_3D,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    
    while size(T0,1) < X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TB = ncon({T0,B},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TB},{[-1,1,2],[1,2,-2,-3,-4]});
        C0 = reshape(M,q*size(T0,1),q*size(T0,1));
        T0 = reshape(TB,[q*size(T0,1),q,q*size(T0,1)]);
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