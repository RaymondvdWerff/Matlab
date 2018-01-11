q = 2;
X = 4;
tol = 10^(-8);

delta = zeros(q,q,q,q);
for i = 1:q
    delta(i,i,i,i) = 1;    
end

spin = zeros(q,q,q,q);
spin(1,1,1,1) = 1;
spin(2,2,2,2) = -1;


T_begin = 1;
T_step = 0.01;
T_end = 3;

emptylist = zeros(floor((T_end-T_begin)/T_step) + 1,1);
T_vals = emptylist;
m_vals1 = emptylist;
m_vals2 = emptylist;
iters = emptylist;
counter = 1;

for Temp = T_begin:T_step:T_end
    
    Q = [exp(1/Temp) exp(-1/Temp); exp(-1/Temp) exp(1/Temp)];
    Qsq = sqrtm(Q);

    T = zeros(q,q,q);
    T(1,1,1) = 1;
    T = ncon({T,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});

    A = ncon({delta,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    B = ncon({spin,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});

    maxiter = 5000;
    minsteps = 4;
    sold = zeros(q,1);

    for i = 1:maxiter

        N = ncon({T,A},{[-1,-3,1],[1,-2,-4,-5]});
        si_N = size(N);
        N_SVD = reshape(N,prod(si_N(1:2)),prod(si_N(3:5)));
        D = min(X,size(N_SVD,1));

        [U,s,~] = svd(N_SVD,'econ');

        U = U(:,1:D);
        U = reshape(U,[si_N(1),si_N(2),D]);
        s = diag(s);
        s = s(1:D);
        snew = s/max(s);

        T = ncon({conj(U),N,U},{[1,2,-1],[1,2,3,4,-3],[3,4,-2]}); 
        T = T+permute(T,[2 1 3]);
        T = T./max(abs(T(:)));

        if numel(sold) == numel(snew)
            diffs = norm(snew-sold);
        else
            diffs = inf;
        end

        if (diffs<tol)  && (i>minsteps) 
            break;
        end
        sold = snew;

    end
    
    M1 = ncon({T,A,conj(T)},{[-1,-4,1],[1,-2,2,-5],[-3,-6,2]});
    M2 = ncon({T,B,conj(T)},{[-1,-4,1],[1,-2,2,-5],[-3,-6,2]});
    si_M = size(M1);
    M_SVD = reshape(M1,prod(si_M(1:3)),prod(si_M(4:6)));

    [U,s,~] = svd(M_SVD,'econ');

    S_k = zeros(size(s,1));
    S_k(1,1) = 1;

    U = reshape(U,[si_M(1),si_M(2),si_M(3),prod(si_M(1:3))]);

    nrm = ncon({M1,U,conj(U),S_k},{[1,2,3,4,5,6],[4,5,6,7],[1,2,3,8],[7,8]});
    m = ncon({M2,U,conj(U),S_k},{[1,2,3,4,5,6],[4,5,6,7],[1,2,3,8],[7,8]});
    m_vals1(counter) = abs(m/nrm);
    
    [V,E] = eigs(reshape(M1,prod(si_M(1:3)),prod(si_M(4:6))),1,'LM');
    V = reshape(V,[si_M(4),si_M(5),si_M(6)]);
    nrm = ncon({V,M1,V},{[1,2,3],[1,2,3,4,5,6],[4,5,6]});
    m = ncon({V,M2,V},{[1,2,3],[1,2,3,4,5,6],[4,5,6]});
    
    m_vals2(counter) = abs(m/nrm);
    T_vals(counter) = Temp;
    iters(counter) = i;
    counter = counter + 1;
    
end

plot(T_vals,m_vals1)
hold on
plot(T_vals,m_vals2,'--')
legend('svd','eig')