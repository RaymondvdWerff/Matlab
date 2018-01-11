%{
D = 2:2:20;
h = 0.8:0.2:1.4;

for j = 1:numel(h)
    EperSite = zeros(numel(D),1);
    for i = 1:numel(D)
        
        tol = 10^(-7);

        maxiter = 2000;
        [HL,HM,HR] = HMPO(h(j));

        eold = 0;
        escenterold = 0;

        N = 2;
        for k = 1:maxiter
            H = ncon({HL,HM,HM,HR},{[-1,-5,1],[-2,-6,1,3],[-4,-8,3,2],[-3,-7,2]});
            N = N + 2;
            [H,ldims,rdims] = lreshape(H,1:4,5:8);
            H = (H+H')/2;

            [Tm,e] = eigs(H,1,'SA');
            escenter = (e - eold)/2;

            Tm = reshape(Tm,ldims);
            [U,s,V] = tensorsvd(Tm,[1,2],[3,4],D(i));

            HL = ncon({HL,U,HM,conj(U)},{[1,4,3],[1,2,-1],[2,5,3,-3],[4,5,-2]});
            HR = ncon({HR,V,HM,conj(V)},{[1,4,3],[1,2,-1],[2,5,-3,3],[4,5,-2]});

            if abs(escenterold-escenter) < tol
                disp('Converged');
                break;
            end

            eold = e;
            escenterold = escenter;
        end
        EperSite(i) = e;
    end
    semilogy(D,EperSite-EperSite(end),'-*')
    hold on
end

xlabel('D')
ylabel('E')
legend(['h = ' num2str(h(1))],['h = ' num2str(h(2))],['h = ' num2str(h(3))],['h = ' num2str(h(4))])
%}

D = [14];
h = 0:0.01:1.4;
Sz = [1,0;0,-1];
Sx = [0,1;1,0];

for i = 1:numel(D)
    Magnetizationx = zeros(numel(h),1);
    Magnetizationz = zeros(numel(h),1);
    for j = 1:numel(h)
        
        tol = 10^(-10);

        maxiter = 2000;
        [HL,HM,HR] = HMPO(h(j));
        
        HL = ncon({HL,diag([1,0]),diag([1,0])},{[1,2,-3],[-1,1],[-2,2]});
        HR = ncon({HR,diag([1,0]),diag([1,0])},{[1,2,-3],[-1,1],[-2,2]});
        
        eold = 0;
        escenterold = 0;

        N = 2;
        for k = 1:maxiter
            H = ncon({HL,HM,HM,HR},{[-1,-5,1],[-2,-6,1,3],[-4,-8,3,2],[-3,-7,2]});
            N = N + 2;
            [H,ldims,rdims] = lreshape(H,1:4,5:8);
            H = (H+H')/2;

            [Tm,e] = eigs(H,1,'SA');
            escenter = (e - eold)/2;

            Tm = reshape(Tm,ldims);
            [U,s,V] = tensorsvd(Tm,[1,2],[3,4],D(i));

            HL = ncon({HL,U,HM,conj(U)},{[1,4,3],[1,2,-1],[2,5,3,-3],[4,5,-2]});
            HR = ncon({HR,V,HM,conj(V)},{[1,4,3],[1,2,-1],[2,5,-3,3],[4,5,-2]});

            if abs(escenterold-escenter) < tol
                disp('Converged');
                break;
            end

            eold = e;
            escenterold = escenter;
        end
        norm = ncon({Tm,conj(Tm)},{[1,2,3,4],[1,2,3,4]});
        mx = ncon({Tm,Sx,Tm},{[2,1,4,5],[1,3],[2,3,4,5]})/norm;
        Magnetizationx(j) = mx;
        mz = ncon({Tm,Sz,Tm},{[2,1,4,5],[1,3],[2,3,4,5]})/norm;
        Magnetizationz(j) = mz;
    end
    plot(h,Magnetizationx)
    hold on
    plot(h,Magnetizationz)
end

xlabel('h')
ylabel('m')
%legend(['D = ' num2str(D(1))],['D = ' num2str(D(2))],['D = ' num2str(D(3))],['D = ' num2str(D(4))])
