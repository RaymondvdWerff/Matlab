function [k,m] = compute_kappa_m(A,Bx,By,C,T)
    
    env_half = ncon({C,C,T,T},{[-1,1],[2,3],[1,-3,2],[3,-4,-2]});
    env = ncon({env_half,env_half},{[1,2,-1,-2],[2,1,-3,-4]});
    kx = ncon({env,Bx},{[1,2,3,4],[1,2,3,4]});
    ky = ncon({env,By},{[1,2,3,4],[1,2,3,4]});
    k1 = ncon({env,A},{[1,2,3,4],[1,2,3,4]});
        
    k2 = ncon({C,C,C,C},{[1,2],[2,3],[3,4],[4,1]});

    CTC = ncon({C,T,C},{[-1,1],[1,-2,2],[2,-3]});
    k3 = ncon({CTC,CTC},{[1,2,3],[1,2,3]});
    
    k = k1*k2/k3^2;
    m = sqrt(kx^2/k1^2+ky^2/k1^2);
    
end

    