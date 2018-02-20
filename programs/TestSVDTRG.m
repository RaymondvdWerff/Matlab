q = 6;
X = 10;
tol = 1e-8;
maxiter = 3;

t_begin = 0.5;
t_step = 0.01;
t_end = 1;

ts = t_begin:t_step:t_end;
emptylist = zeros(1,numel(ts));
XX = emptylist;
iters = emptylist;
tictocs = emptylist;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end

for t = 1:numel(ts)
    
    temp = ts(t);
    disp(['temp = ' num2str(temp)]);
    
    tic
    Qsq = sqrtm(Q_clock(q,temp,0));
    A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    
    for iter = 1:maxiter
        
        [S1,s1,S3] = tensorsvd(A,[1,3],[2,4],X,'m');
        s1 = s1/max(s1(:));
        
        [S2,s2,S4] = tensorsvd(A,[3,2],[4,1],X,'m');
        s2 = s2/max(s2(:));
        
        wi = eye(size(S1,1));
        wj = eye(size(S1,1));
        wk = eye(size(S1,1));
        wl = eye(size(S1,1));
        
        for i = 1:100
            
            M12 = ncon({sqrtm(wi),S2,S1,sqrtm(wk)},{[-1,1],[1,2,-2],[2,3,-4],[-3,3]}); 
            [U,D12,V] = tensorsvd(M12,[1,2],[3,4],min(size(S1,1),X));
            
            S1 = ncon({inv(sqrtm(wk)),V,sqrtm(D12)},{[-2,1],[1,-3,2],[2,-1]});
            S1 = S1/max(abs(S1(:)));
            S2 = ncon({inv(sqrtm(wi)),U,sqrtm(D12)},{[-1,1],[1,-3,2],[2,-2]});
            S2 = S2/max(abs(S2(:)));
            
            M23 = ncon({sqrtm(wl),S3,S2,sqrtm(wj)},{[-1,1],[1,2,-2],[2,3,-4],[-3,3]}); 
            [U,D23,V] = tensorsvd(M23,[1,2],[3,4],min(size(S2,1),X));
            
            S2 = ncon({inv(sqrtm(wl)),V,sqrtm(D23)},{[-2,1],[1,-3,2],[2,-1]});
            S2 = S2/max(abs(S2(:)));
            S3 = ncon({inv(sqrtm(wl)),U,sqrtm(D23)},{[-1,1],[1,-3,2],[2,-2]});
            S3 = S3/max(abs(S3(:)));
            
            M34 = ncon({sqrtm(wk),S4,S3,sqrtm(wi)},{[-1,1],[1,2,-2],[2,3,-4],[-3,3]}); 
            [U,D34,V] = tensorsvd(M34,[1,2],[3,4],min(size(S3,1),X));
            
            S3 = ncon({inv(sqrtm(wi)),V,sqrtm(D34)},{[-2,1],[1,-3,2],[2,-1]});
            S3 = S3/max(abs(S3(:)));
            S4 = ncon({inv(sqrtm(wk)),U,sqrtm(D34)},{[-1,1],[1,-3,2],[2,-2]});
            S4 = S4/max(abs(S4(:)));
            
            M41 = ncon({sqrtm(wj),S1,S4,sqrtm(wl)},{[-1,1],[1,2,-2],[2,3,-4],[-3,3]}); 
            [U,D41,V] = tensorsvd(M41,[1,2],[3,4],min(size(S4,1),X));
            
            S4 = ncon({inv(sqrtm(wl)),V,sqrtm(D41)},{[-2,1],[1,-3,2],[2,-1]});
            S4 = S4/max(abs(S4(:)));
            S1 = ncon({inv(sqrtm(wj)),U,sqrtm(D41)},{[-1,1],[1,-3,2],[2,-2]});
            S1 = S1/max(abs(S1(:)));
            
            wi = D23;
            wj = D12;
            wk = D41;
            wl = D34;
        end
        
        A = ncon({S1,S2,S3,S4},{[1,4,-4],[2,1,-1],[3,2,-3],[4,3,-2]});
        A = A/max(A(:));
        if iter > 10
            delta = sum(sum(abs(s1-s1old)));
            if delta < tol
                break;
            end
        end  
        s1old = s1;
     
    end
    
    if (iter == maxiter) && (maxiter > 10)
        disp(['SVDTRG not converged at T = ' num2str(temp)]);
    end 
    
    X1 = ncon({A},{[1,1,2,2]});
    X2 = ncon({A,A},{[1,2,3,4],[2,1,4,3]});

    XX(t) = X1^2/X2;
    tictocs(t) = toc;
    iters(t) = iter;  
end

plot(ts,XX,'o-');