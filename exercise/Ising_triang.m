q = 2;

delta = zeros(q,q,q,q,q,q);
for i = 1:q
    delta(i,i,i,i,i,i) = 1;    
end

spin = zeros(q,q,q,q,q,q);
spin(1,1,1,1,1,1) = 1;
spin(q,q,q,q,q,q) = -1;

tol = 10^(-6);
b_begin = 0.2; 
b_end = 1.2;
b_step = 0.001;

X = [10];
emptylist = zeros(floor((b_end-b_begin)/b_step) + 1,1);

for x = 1:numel(X)
    
    m_vals = emptylist;
    T_vals = emptylist;
    exact_vals = emptylist;
    iters = emptylist;

    counter = 1;

    for b = b_begin:b_step:b_end

        %Q = ones(q);
        %for i = 1:q
        %    Q(i,i) = exp(b);
        %end
        Q = [exp(b) exp(-b); exp(-b) exp(b)];
        Qsq = sqrtm(Q);

        A = ncon({delta,Qsq,Qsq,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4,5,6],[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-6,6]});
        B = ncon({spin,Qsq,Qsq,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4,5,6],[-1,1],[-2,2],[-3,3],[-4,4],[-5,5],[-6,6]});

        [C1,C2,T1,T2,iter] = renorm_Ising_triang(Qsq,X(x),A,tol);
        m = abs(collapse_triang(C1,C2,T1,T2,B)/collapse_triang(C1,C2,T1,T2,A));
        %m = (q*n-1)/(q-1);

        m_vals(counter) = m;
        T_vals(counter) = 1/b;
        iters(counter) = iter;

        counter = counter + 1;
    end
    plot(T_vals,m_vals)
    hold on
end
title('Magnetization as a function of temperature for the triangular Ising model')
xlabel('Temperature')
ylabel('Magnetization')
%legend(['X = ' num2str(X(1))],['X = ' num2str(X(2))],['X = ' num2str(X(3))],['X = ' num2str(X(4))])
