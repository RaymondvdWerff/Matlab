function [Q,L] = ql(A)

    [Qi,Ri] = qr(fliplr(A),0);

    L = fliplr(flipud(Ri));
    Q = fliplr(Qi);
    
end