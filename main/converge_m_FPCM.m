function [m,sv,iters,tictocs,imarkers,tmarkers] = converge_m_FPCM(Q,q,X,tol,maxiter,temp,tols)
    
    m = zeros(1,maxiter);
    sv = zeros(1,maxiter);
    iters = 1:maxiter;
    p = 12;
    tol2 = 1e-3;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    %spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    
    Qsq = sqrtm(Q(q,temp,0));Qsq = real(Qsq);
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    %B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
    [~,T] = beginmatrices(Qsq,A,X,1);
    tictocs = [0];
    i = 1;
      
    for iter = 1:maxiter
        tic;
        [Tl,C] = LeftOrthonormalize(T,min(1e-6,tol),maxiter/10,temp,p,tol2);
        [~,s,~] = svd(C);s = s/max(s(:));
        
        opts.v0 = reshape(T,q*X^2,1);
        opts.tol = tol2;
        opts.p = p;
        [T,~] = eigs(@(x)mult2(Tl,A,x),q*X^2,1,'LM',opts);
        T = reshape(T,[X,q,X]);
        T = T + permute(T,[3,2,1]);
        
        tictocs(iter) = tictocs(end) + toc;
        
        env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
        env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
        Z = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        mx = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]})/Z;
        my = ncon({env,By},{[1,2,3,4],[1,2,3,4]})/Z;
        
        m(iter) = sqrt(mx^2+my^2);
        
        if iter > 1
            sv(iter) = sum(sum(abs(s-sold)));
            if i < numel(tols)+1
                if sum(sum(abs(s-sold))) < tols(i)
                    imarkers(i) = iter;
                    tmarkers(i) = tictocs(iter);
                    i = i+1;
                end
            end
%             if abs(m(iter)-m(iter-1)) < tol
%                 break;
%             end
        end
        sold = s;
    end
    
%     if iter == maxiter
%         disp('converge_m_FPCM not converged');
%     end 
%     iters = 1:iter;
end

function [Tl,C1] = LeftOrthonormalize(T,tol,maxiter,temp,p,tol2)

    X = size(T,1);
    q = size(T,2);
    
    for iter = 1:maxiter
        
        if iter == 1
            C2 = ncon({T,T},{[1,2,-2],[1,2,-1]});
            opts.v0 = reshape(C2,X^2,1);
            opts.tol = tol2;
            opts.p = p;
            [C2,~] = eigs(@(x)mult1(T,T,x),X^2,1,'LM',opts);
            C2 = reshape(C2,X,X);

            [U,s2,~] = svd(C2);
            C = U*sqrt(s2)*U';
            C = C./max(abs(C(:)));       
        else
            opts.v0 = reshape(C1,X^2,1);
            opts.tol = tol2;
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
    y = ncon({T,M1,M2,M1},{[4,2,1],[1,3,-3],[3,5,2,-2],[4,5,-1]});
    y = reshape(y,prod(si),1);
end
