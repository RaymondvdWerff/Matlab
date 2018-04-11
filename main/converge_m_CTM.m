function [m,sv,iters,tictocs,imarkers,tmarkers] = converge_m_CTM(Q,q,X,tol,maxiter,temp,tols)
    
    m = zeros(1,maxiter);
    sv = zeros(1,maxiter);
    iters = 1:maxiter;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    %spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    
    Qsq = sqrtm(Q(q,temp,0));Qsq = real(Qsq);
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    %B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    
    [C,T] = beginmatrices(Qsq,A,X,1);
    
    i = 1;
    tictocs = [0];
       
    for iter = 1:maxiter
        
        tic
        
        N = ncon({T,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({C,T,N},{[1,2],[1,3,-1],[2,3,-2,-3,-4]});
        
        [U,s,~] = tensorsvd(M,[1,2],[3,4],X);s = s/max(s(:));
    
        C = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C = (C + permute(C,[2,1]))/max(C(:));
        
        T = ncon({U,N,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T = (T + permute(T,[3,2,1]))/max(T(:));
        
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
%         disp('converge_m_CTM not converged');
%     end 
%     iters = 1:iter;
end
