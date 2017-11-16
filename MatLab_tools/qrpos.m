function [Q,R] = qrpos(A)
    [Q,R] = qr(A,0);
    signR = diag(sign(diag(R)));
    Q = Q*signR;
    R = signR*R;
end
    