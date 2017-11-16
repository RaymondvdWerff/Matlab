function E_tot = compute_E(MPS,H)
    
    E_tot = 0;
    
    for n = 1:numel(MPS)-1
        E = ncon({MPS{n},MPS{n+1},H,conj(MPS{n}),conj(MPS{n+1})},{[7,3,1],[1,4,8],[3,4,5,6],[7,5,2],[2,6,8]});
        norm = tcontract(MPS{n},conj(MPS{n}),[1,2,3],[1,2,3]);
        E_tot = E_tot + E/norm;
        
        [U,s,V] = tensorsvd(MPS{n},[1,2],3,inf,'r');
        MPS{n} = U;
        MPS{n+1} = tcontract(permute(V,[2,1]),MPS{n+1},2,1);
        
    end
end
        