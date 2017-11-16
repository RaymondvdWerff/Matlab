function H = get_H_trans_ising(lambda)
    % Returns the 2-site Hamiltonian of the ferromagnetic transverse field Ising
    % model in 2D
    % The site-term is included in the 2-site operator in a symmetric way
    sigmaz = [[1 0];[0 -1]];
    sigmax = [[0 1];[1 0]];
    id = eye(2);
    H = -kron(sigmaz,sigmaz)-lambda/4*(kron(sigmax,id)+kron(id,sigmax));
    H = reshape(H,2,2,2,2);
end