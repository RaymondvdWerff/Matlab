function E = Energy(H,C,T,A,D)
    
    s = size(T);
    T = reshape(T,[s(1),D,D,s(3)]);
    
    half = ncon({C,C,T,T,T,A,conj(A)},{[1,3],[2,10],[-1,5,7,3],[1,6,4,2],[10,8,9,-4],[4,5,8,-2,-5],[6,7,9,-3,-6]});
    rho = ncon({half,half},{[1,2,3,4,-1,-2],[1,2,3,4,-3,-4]});
    norm = ncon({rho},{[1,1,2,2]});
    E = ncon({rho,H},{[1,2,3,4],[1,2,3,4]});
    E = E/norm;
    
end
