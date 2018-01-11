function [A_new,B_new] = contract(A_old,B_old,M,D,lr)

    contr = ncon({A_old,B_old,M},{[-1,2,1],[1,3,-4],[2,3,-2,-3]});
    [U,s,V] = tensorsvd(contr,[1,2],[3,4],D,lr);
    
    A_new = U;
    B_new = permute(V,[3,1,2]);
    
    norm = trace(s*s);
    
    if lr == 'l'
        A_new = A_new./sqrt(norm);
    else    
        B_new = B_new./sqrt(norm);
    end
end       