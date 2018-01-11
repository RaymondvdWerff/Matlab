function Eres = ex_quantum_ising(lambda)
    % Perform a (brute-force) energy minimzation for the
    % transverse-field 2D Ising model
    
    if nargin<1 
        %lambda=3.04438; critical point
        lambda=3;
    end
    
    % get Hamiltonian
    H = get_H_trans_ising(lambda);
    
    % initialize random symmetric tensor using 12 random parameters
    disp('Test: random tensor');
    rand('seed',2);
    c = rand(12,1)-0.5;
    A = get_symtensor(c);
    chi = 8;
    tol=1e-6;
    % compute energy
    E=doctmq(A, H, chi, tol, 1);
    fprintf('Energy of random initial tensor: %g',E);
    
    % do a brute-force optimization to minimize energy
    function E=get_E(c)
        A = get_symtensor(c);
        E=doctmq(A,H,chi,tol,0);
    end
    
    x1 = ones(1,12);
    opts.TolFun=1e-8;
    opts.MaxFunEvals=10000;
    [cres,Eres] = fmincon(@get_E,c,[],[],[],[],-x1,x1,[],opts);
    disp('Results:');
    fprintf('lambda=%g, Energy (D=2):%g\n', lambda, Eres);
end
