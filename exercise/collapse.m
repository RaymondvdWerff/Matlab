function B = collapse(C,T,A)
    CTC = ncon({C,T,C},{[-1,1],[1,-2,2],[2,-3]});
    TAT = ncon({T,A,T},{[-1,1,-4],[1,-2,2,-5],[-3,2,-6]});
    upper = ncon({CTC,TAT},{[1,2,3],[1,2,3,-1,-2,-3]});
    B = ncon({CTC,upper},{[1,2,3],[1,2,3]});
end