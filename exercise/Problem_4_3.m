l_begin = 1; 
l_end = 4;
l_step = 0.1;

emptylist = zeros(floor((l_end-l_begin)/l_step) + 1,1);
ls = emptylist;
Es = emptylist;
tol = 10^(-7);
X = 2;
counter = 1;

for l = l_begin:l_step:l_end
    A = get_symtensor([1 0.5 0.3 0.5 0.1 0.1 1 0.2 0.3 0.5 0.1 1]);
    D = size(A,1);
    Qsq = rand(D^2);
    H = get_H_trans_ising(l);

    AA = ncon({A,conj(A)},{[-1,-2,-3,-4,1],[-5,-6,-7,-8,1]});
    AA = reshape(AA,[D^2,D^2,D^2,D^2]);

    [C,T,i] = renorm_Ising2(Qsq,X,AA,tol);

    Es(counter) = Energy(H,C,T,A,D);
    %Es(counter) = doctmq(A, H, X, tol, 0);
    ls(counter) = l;
    counter = counter + 1;
end

plot(ls,Es,'*')