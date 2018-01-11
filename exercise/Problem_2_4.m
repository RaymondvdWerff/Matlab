L = 16;
J = 1;
tau = 0.1;

Sz = diag([0.5 -0.5]);
Sp = [0 1;0 0];
Sm = Sp';
    
H2S = J*(kron(Sz,Sz)+0.5*(kron(Sp,Sm)+kron(Sm,Sp)));

T1 = reshape(expm(-tau*H2S),2,2,2,2);
H = reshape(H2S,2,2,2,2);

Energys = [];
Ds = [];

for D = 2:2:20
    state = {randn(1,2,D)};
    for i = 2:L-1
        state{i} = randn(D,2,D);
    end
    state{L} = randn(D,2,1);

    Energy = [];
    Iteration = [];
    E_old = 1;

    for i = 1:200
        for n = 1:(L-1)
            [state{n},state{n+1}] = contract(state{n},state{n+1},T1,D,'r');
        end

        for n = 1:(L-1)
            [state{L-n},state{L+1-n}] = contract(state{L-n},state{L+1-n},T1,D,'l');
        end

        if mod(i,50) == 0
            E = compute_E(state,H);
            Energy = [Energy,E];
            Iteration = [Iteration,i];
        end

        if abs(E_old-E) < 10^(-5)
            disp('Converged')
            break

        E_old = E;

        end

    end

    Energys = [Energys,Energy(end)];
    Ds = [Ds,D];
    
end

plot(1./Ds,Energys,'*')
    
