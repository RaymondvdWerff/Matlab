function H = Hamiltonian(N)
    J = 1;
    
    Sz = diag([0.5 -0.5]);
    Sp = [0 1;0 0];
    Sm = Sp';
    
    H2S = J*(kron(Sz,Sz)+0.5*(kron(Sp,Sm)+kron(Sm,Sp)));
    H2S = sparse(H2S);
    H = H2S;
    
    for i = 1:(N-2)
        H = kron(H,speye(2))+kron(speye(2^i),H2S);
    end
end

