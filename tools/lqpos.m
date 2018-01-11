function [L,Q] = lqpos(A)

    [Qi, Li]=ql(flipud(A).');
    L = fliplr(flipud(Li.'));
    Q = flipud(Qi.');
    signL = diag(sign(diag(L)));
    L = L*signL;
    Q = signL*Q;

end
