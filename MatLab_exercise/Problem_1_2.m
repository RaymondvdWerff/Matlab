L = 8;

state0 = rand(2^L,1);
state = reshape(state0,2,2^(L-1));
MPS = cell(L,1);

for i = 1:L-1
    
    [U,s,V] = svd(state);
    s1 = size(s,1);
    s2 = size(s,2);
    
    if i == 1
        MPS{i} = U;
    else
        MPS{i} = reshape(U,s1_old,2,s1);
    end
    
    state = s*V';
    state = reshape(state,s1*2,s2/2);
    s1_old = s1;
    
end

MPS{L} = reshape(state,s1,s2);

state1 = MPS{1};

for i = 2:L-1
    state1 = tcontract(state1,MPS{i},2,1);
    state1 = reshape(state1,2^i,size(MPS{i},3));
end

state1 = tcontract(state1,MPS{L},2,1);
state1 = reshape(state1,2^L,1);

E0 = state0'*Hamiltonian(L)*state0/(state0'*state0)
E1 = state1'*Hamiltonian(L)*state1/(state1'*state1)





    