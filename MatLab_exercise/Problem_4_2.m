q = 2;

delta = zeros(q,q,q,q);
for i = 1:q
    delta(i,i,i,i) = 1;    
end

spin = zeros(q,q,q,q);
spin(1,1,1,1) = 1;

tol = 10^(-6);
b_begin = 0.8; 
b_end = 1.2;
b_step = 0.001;

X = [2];
emptylist = zeros(floor((b_end-b_begin)/b_step) + 1,1);

for x = 1:numel(X)
    
    m_vals = emptylist;
    T_vals = emptylist;
    exact_vals = emptylist;
    iters = emptylist;

    counter = 1;

    for b = b_begin:b_step:b_end

        Q = ones(q);
        for i = 1:q
            Q(i,i) = exp(b);
        end
        Qsq = sqrtm(Q);

        A = ncon({delta,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        B = ncon({spin,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});

        [C,T,iter] = renorm_Ising2(Qsq,X(x),A,tol);
        n = abs(collapse(C,T,B)/collapse(C,T,A));
        m = (q*n-1)/(q-1);

        m_vals(counter) = m;
        T_vals(counter) = 1/b;
        iters(counter) = iter;

        counter = counter + 1;
    end
    plot(T_vals,m_vals)
    hold on
end
title(['Magnetization as a function of temperature for the ' num2str(q) '-state Potts model'])
xlabel('Temperature')
ylabel('Magnetization')
%legend(['X = ' num2str(X(1))],['X = ' num2str(X(2))],['X = ' num2str(X(3))],['X = ' num2str(X(4))])
