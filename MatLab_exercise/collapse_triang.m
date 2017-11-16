function B = collapse_triang(C1,C2,T1,T2,A)
    half = ncon({C1,T1,C2,T2},{[-1,1],[1,2,-2,-3],[2,3,-4],[3,-7,-6,-5]});
    B = ncon({half,half,A},{[4,5,1,2,3,6,7],[7,6,8,9,10,5,4],[1,2,3,8,9,10]});
end