function TestEigs
    d = 2;
    D = 50;
    A = rand(d,D,D);
    A = A + permute(A,[1,3,2]);
    tic
    opts.v0 = ones(D^2,1);
    [V1,D1] = eigs(@(x)mult(x,A,D),D^2,1,'LM',opts);
    D1
    toc
    tic
    M = ncon({A,conj(A)},{[1,-1,-3],[1,-2,-4]});
    M = reshape(M,D^2,D^2);
    [V2,D2] = eigs(M,1,'LM');
    D2
    toc
end

function y = mult(x,A,D)
    x = reshape(x,D,D);
    y = ncon({A,conj(A),x},{[2,-1,1],[2,-2,3],[1,3]});
    y = reshape(y,D^2,1);
end