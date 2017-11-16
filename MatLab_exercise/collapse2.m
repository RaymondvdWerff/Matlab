function B = collapse2(C,T)
    CTC = ncon({C,T,C},{[-1,1],[1,-2,2],[2,-3]});
    B = ncon({CTC,CTC},{[1,2,3],[1,2,3]});
end