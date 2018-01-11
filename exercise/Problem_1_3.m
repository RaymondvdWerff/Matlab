L = 10;

state0 = rand(2^L,1);
E0 = state0'*Hamiltonian(L)*state0/(state0'*state0);

Errors = [];
Ds = [];

for D = 2:100
    
    state = reshape(state0,2,2^(L-1));
    MPS = cell(L,1);

    for i = 1:L-1

        [U,s,V] = svd(state);
        X = min([D size(s,1) size(s,2)]);

        if i == 1
            MPS{i} = U(:,1:X);
        else
            MPS{i} = reshape(U(:,1:X),s1_old,2,X);
        end

        state = s(1:X,1:X)*V(:,1:X)';
        a = size(state,1);
        b = size(state,2);
        state = reshape(state,a*2,b/2);
        s1_old = a;

    end

    MPS{L} = reshape(state,X,2);

    state1 = MPS{1};

    for i = 2:L-1
        state1 = tcontract(state1,MPS{i},2,1);
        state1 = reshape(state1,2^i,size(MPS{i},3));
    end

    state1 = tcontract(state1,MPS{L},2,1);
    state1 = reshape(state1,2^L,1);

    E1 = state1'*Hamiltonian(L)*state1/(state1'*state1);
    Errors = [Errors,abs(E0-E1)];
    Ds = [Ds,D];
    
end

semilogy(Ds,Errors)
xlabel('Bond dimension')
ylabel('Relative error')
  
