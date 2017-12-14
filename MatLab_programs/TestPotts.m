%Testing the Potts model algorithms.
%{
%Plot magnetization vs temperature.
q = 2;
X = 20;
tol1 = 1e-4;
tol2 = 1e-4;
maxiter = 10000;

t_begin = 1.0;
t_step = 1e-3;
t_end = 1.2;

ts = t_begin:t_step:t_end;
m = m_exact(ts);

[m1,iters1,tictocs1] = compute_m(q,X,tol1,maxiter,ts,@Potts_CTM);
[m2,iters2,tictocs2] = compute_m(q,X,tol2,maxiter,ts,@Potts_FPCM);

plot(ts,tictocs1,'+-');hold on;
plot(ts,tictocs2,'x-');hold on;
%plot(ts,m,'o-');

title(['Magnetization as a function of temperature for the ' num2str(q) '-state Potts model']);
xlabel('temperature');
ylabel('magnetization');
%ylabel('time (s)');
%axis([ts(1),ts(end),-0.1,1]);
%legend(['X = ' num2str(X(1))],['X = ' num2str(X(2))],['X = ' num2str(X(3))],['X = ' num2str(X(4))],['X = ' num2str(X(5))],'Exact')
%legend(['tol = ' num2str(tols(1))],['tol = ' num2str(tols(2))],['tol = ' num2str(tols(3))],['tol = ' num2str(tols(4))],['tol = ' num2str(tols(5))])
legend('CTM','FPCM');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot correlation length vs temperature.
q = 2;
X = 20;
tols = [1e-3];
maxiter = 10000;

t_begin = 1.13;
t_step = 1e-6;
t_end = 1.14;

ts = t_begin:t_step:t_end;

for tol = 1:numel(tols)
    [corr,iters,tictocs] = compute_corr(q,X,tols(tol),maxiter,ts,@Potts_FPCM);
    plot(ts,corr,'o-');hold on;
end

title(['Correlation length as a function of temperature for the ' num2str(q) '-state Potts model']);
xlabel('temperature');
ylabel('correlation length');
legend(['tol = ' num2str(tols(1))],['tol = ' num2str(tols(2))],['tol = ' num2str(tols(3))],['tol = ' num2str(tols(4))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test FPCM for different Krylov subspace dimensions.
%{
q = 2;
X = 4;
tol = 1e-4;
maxiter = 5000;

t_begin = 1.1;
t_step = 1e-3;
t_end = 1.2;

ts = t_begin:t_step:t_end;
K = 4:4:q*X^2;

for k = 1:numel(K)
    K(k)
    tic
    [m,iters,tictocs] = compute_m(q,X,tol,maxiter,K(k),ts,@Potts_FPCM);
    toc
    plot(ts,tictocs,'o-');hold on;
end

title('Magnetization as a function of temperature for different Krylov subspace dimensions');
xlabel('temperature');
ylabel('magnetization');
legend(['k = ' num2str(K(1))],['k = ' num2str(K(2))],['k = ' num2str(K(3))],['k = ' num2str(K(4))],['k = ' num2str(K(5))],['k = ' num2str(K(6))],['k = ' num2str(K(7))],['k = ' num2str(K(8))])
%}   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine magnetization as a function of iteration number/time
q = 2;
X = 20;
%tc = 1/log(1+sqrt(2));
%t = tc + 1e-3;
temp = 1.135;
maxiter1 = 5000;
maxiter2 = 1000;
tol = 1e-8;
tols = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];

[m1,iters1,tictocs1,imarkers1,tmarkers1] = converge_m_CTM(q,X,tol,maxiter1,temp,tols);
[m2,iters2,tictocs2,imarkers2,tmarkers2] = converge_m_FPCM(q,X,tol,maxiter2,temp,tols);
%{
%As a function of iteration number.
plot(iters1,m1,'b');hold on;
plot(iters1(end),m1(end),'*b');hold on;
plot(iters2,m2,'r');hold on;
plot(iters2(end),m2(end),'*r');hold on;
for i = 1:numel(imarkers1)
    plot([imarkers1(i),imarkers1(i)],[min(min(m1),min(m2)),max(max(m1),max(m2))],'b--');hold on;
end
for i = 1:numel(imarkers2)
    plot([imarkers2(i),imarkers2(i)],[min(min(m1),min(m2)),max(max(m1),max(m2))],'r--');hold on;
end
xlabel('iteration number');
ylabel('Magnetization');
title(['Magnetization as a function of iteration number for the ' num2str(q) '-state Potts model for X = ' num2str(X) ' at T = ' num2str(temp)]);
legend('CTM','','FPCM','');
%}

%As a function of time.
plot(tictocs1,m1,'b');hold on;
plot(tictocs1(end),m1(end),'*b');hold on;
plot(tictocs2,m2,'r');hold on;
plot(tictocs2(end),m2(end),'*r');hold on;
for i = 1:numel(tmarkers1)
    plot([tmarkers1(i),tmarkers1(i)],[min(min(m1),min(m2)),max(max(m1),max(m2))],'b--');hold on;
end
for i = 1:numel(tmarkers2)
    plot([tmarkers2(i),tmarkers2(i)],[min(min(m1),min(m2)),max(max(m1),max(m2))],'r--');hold on;
end
xlabel('time (s)');
ylabel('Magnetization');
title(['Magnetization as a function of time for the ' num2str(q) '-state Potts model for X = ' num2str(X) ' at T = ' num2str(temp)]);
legend('CTM','','FPCM','');

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine the critical temperature for different bond dimensions.
format long

q = 2;
X = [20];
tol = 1e-3;
maxiter = 5000;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
tcs = zeros(numel(X),1);

for x = 1:numel(X)
    disp(['X = ' num2str(X(x))]);
    
    t_begin = 1.0;
    t_step = 0.1;
    t_end = 1.5;

    tic
    for l = 1:6
        
        ts = t_begin:t_step:t_end;
        emptylist = zeros(numel(ts),1);
        m = emptylist;
        corr = emptylist;
        difference = 0;
        i = 1;
        
        for t = t_begin:t_step:t_end
            
            Qsq = sqrtm(ones(q)+(exp(1/t)-1)*eye(q));
            A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
            B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
            [C,T] = beginmatrices(Qsq,A,X(x));
            
            [C,T,~] = Potts_FPCM(A,C,T,q,X(x),tol,maxiter,t);

            %n = collapse(C,T,B)/collapse(C,T,A);
            %m(i) = (q*n-1)/(q-1);

            %if i > 1
            %    if (abs(m(i)-m(i-1)))>difference
            %        difference = (abs(m(i)-m(i-1)));
            %        tc = t;
            %        index = i;
            %    
            %    else
            %        break;
            %    end
            %end
            
            corr(i) = corrlen(T);
            if i > 1
                if corr(i) > corr(i-1)
                    tc = t;
                    index = i;
                    
                else
                    break;
                end
            end
            
            i = i + 1; 
        end
        t_begin = ts(index)-t_step;
        t_end = ts(index)+t_step;
        t_step = t_step/10;
    end
    toc
    tc
    tcs(x) = tc;
end
%}
%{
%Alg tol=10^(-4) t_step=10^(-5)
%X1 = [2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75];
%tcs1 = [1.16333,1.15806,1.14059,1.13989,1.13720,1.13686,1.13703,1.13617,1.13550,1.13503,1.13492,1.13474,1.13471,1.13468,1.13466,1.13465,1.13464,1.13463,1.13462,1.13462,1.13462,1.13462];
%comptime1 = [236,272,648,270,676,655,729,334,897,790,960,1392,1987,2594,3483,4499,6846,9788,9181,9987,13168,13569];

%Alg tol=10^(-4) t_step=10^(-6)
%X1 = [20,25,30,35,40,45,50,55,60,65,70,75];
%tcs1 = [1.134902,1.134738,1.134720,1.134672,1.134656,1.134647,1.134633,1.134630,1.134624,1.134619,1.134613,1.134613];
%comptime1 = [1336,2131,3082,3956,5077,7203,10194,12856,15983,18138,23728,27544];

%CTM2 tol=10^(-4) t_step=10^(-5)
%X1 = [2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75];
%tcs1 = [1.16240,1.15568,1.14090,1.13927,1.13727,1.13673,1.13660,1.13649,1.13571,1.13505,1.13552,1.13549,1.13550,1.13547,1.13547,1.13547,1.13547,1.13547,1.13546,1.13547,1.13547,1.13546];
%comptime1 = [28,25,55,45,65,71,69,72,94,114,123,138,153,174,196,227,266,308,352,415,475,568];

%CTM2 tol=10^(-8) t_step=10^(-5)
%X1 = [2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75];
%tcs1 = [1.16320,1.15800,1.14049,1.13988,1.13713,1.13614,1.13686,1.13609,1.13549,1.13505,1.13492,1.13485,1.13482,1.13480,1.13479,1.13478,1.13478,1.13478,1.13478,1.13478,1.13478,1.13477];    
%comptime1 = [420,383,429,379,463,444,471,472,468,519,579,618,688,762,836,960,1092,1235,1366,1547,1718,1944];

%CTM2 tol=10^(-8) t_step=10^(-6)
%X1 = [20,25,30,35,40,45,50,55,60,65,70,75];
%tcs1 = [1.134912,1.134847,1.134815,1.134797,1.134782,1.134779,1.134776,1.134773,1.134772,1.134774,1.134772,1.134770];    
%comptime1 = [925,916,1008,1125,1229,1395,1572,1784,1951,2242,2452,2754];

%FPCM tol=10^(-4) t_step=10^(-6)
X1 = [20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100];
tcs1 = [1.134842,1.134740,1.134716,1.134670,1.134656,1.134648,1.134633,1.134622,1.134623,1.134611,1.134614,1.134612,1.134607,1.134609,1.134607,1.134605,1.134604];
coomptime1 = [209,450,218,748,628,908,1034,712,1001,843,1595,1624,1774,2512,2741,2804,4600];

begin = 1;
X2 = reshape(X1(begin:end),numel(X1)-begin+1,1);
tcs2 = reshape(tcs1(begin:end),numel(tcs1)-begin+1,1);

X3 = 20:0.01:100;
X4 = 10:0.01:100;

myfittype = fittype('a + b*x^c');
myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1.1346,0.1,-2]);
myfit = fit(X2,tcs2,myfittype,myfittopts);
coeff = coeffvalues(myfit);Tc=coeff(1);b=coeff(2);c=coeff(3);
err_coeff = confint(myfit);err_Tc = err_coeff(1);

subplot(1,2,1);
plot(X3,Tc+b*X3.^c,'r');hold on;
plot(X1,tcs1,'ob');
plot(X2,tcs2,'ob','MarkerFaceColor','b');
xlabel('X');
ylabel('Tc''');
legend('fit','data','data used for fit');
title('Critical temperature as a function of bond dimension for the Alg-method');

subplot(1,2,2);
loglog(1./X4,b*X4.^(c),'r');hold on;
loglog(1./X1,abs(tcs1-Tc),'ob');
loglog(1./X2,abs(tcs2-Tc),'ob','MarkerFaceColor','b');
xlabel('1/X');
ylabel('Tc''-Tc*');
legend('fit','data','data used for fit');
title('Difference between Tc''(X) and Tc* determined from fit as a function of 1/X');

disp('Critical temperature from fit:');
Tc
error = abs(Tc-err_Tc)
disp('Exact critical temperature:');
Tc_exact = 1/log(sqrt(2)+1)
%}
%{
%Plot of the computation time as function of bond dimension.
begin = 1;
X2 = reshape(X1(begin:end),numel(X1)-begin+1,1);
comptime2 = reshape(comptime1(begin:end),numel(comptime1)-begin+1,1);

myfittype2 = fittype('a+b*x^c');
%myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1.1346,0.1,-2]);
myfit2 = fit(X2,comptime2,myfittype2)

plot(myfit2,X2,comptime2);
xlabel('X');
ylabel('time (s)');
title('Computation time of the critical temperature as a function of bond dimension');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Plot free energy vs temperature.
q = 2;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
spin1 = zeros(q,q,q,q);spin1(1,1,1,1)=1;

t_begin = 0.5;
t_step = 0.05;
t_end = 3;

ts = t_begin:t_step:t_end;
emptylist = zeros(numel(ts),1);

f3 = emptylist;

X = [4];

tol1 = 10^(-6);
tol2 = 10^(-4);

for x = 1:numel(X) 
    disp(['X = ' num2str(X(x))]);
    
    tic
    f1 = emptylist;
    f2 = emptylist;
    f3 = emptylist;
    iters1 = emptylist;
    iters2 = emptylist;
    i = 1;
    
    T1 = rand(X(x),q,X(x));T1 = T1 + permute(T1,[3,2,1]);
    C1 = rand(X(x),X(x));C1 = C1 + permute(C1,[2,1]);
    T2 = T1;
    C2 = C1;
    
    for t = t_begin:t_step:t_end
        
        Qsq = sqrtm(ones(q)+(exp(1/t)-1)*eye(q));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        B = ncon({spin1,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
        [T1,C1,iter1] = Potts_CTM(T1,C1,B,A,q,X(x),tol1,t);
        [T2,C2,iter2] = Potts_Alg(T2,C2,B,A,q,X(x),tol2,t);
        
        k1 = collapse(C1,T1,A)*ncon({C1,C1,C1,C1},{[1,2],[2,3],[3,4],[4,1]})/collapse2(C1,T1)^2;
        k2 = collapse(C2,T2,A)*ncon({C2,C2,C2,C2},{[1,2],[2,3],[3,4],[4,1]})/collapse2(C2,T2)^2;
        f1(i) = -t*log(k1);
        f2(i) = -t*log(k2);
        
        if x==1;f3(i)=f_exact(t);end
        iters1(i) = iter1;
        iters2(i) = iter2;
        i = i + 1; 
    
    end
        
    toc
    
    plot(ts,f1,'*-');
    hold on
    plot(ts,f2,'o-');
    hold on
end

plot(ts,f3,'--');
title(['Free energy per site as a function of temperature for the ' num2str(q) '-state Potts model'])
xlabel('temperature')
ylabel('free energy')
%legend(['X = ' num2str(X(1))],['X = ' num2str(X(2))],['X = ' num2str(X(3))],['X = ' num2str(X(4))])
%legend(['tol = ' num2str(tols(1))],['tol = ' num2str(tols(2))],['tol = ' num2str(tols(3))],['tol = ' num2str(tols(4))],['tol = ' num2str(tols(5))])
legend('CTM','Alg','Exact');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Test LeftOrhonormalize6
q = 2;
temp = 0.4;
X = 4;
tol = 1e-6;

Qsq = sqrtm(ones(q)+(exp(1/temp)-1)*eye(q));
A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
[~,T] = beginmatrices(Qsq,A,X);

[Tl,C] = LeftOrthonormalize6(T,tol,temp);

T1 = ncon({C,T},{[-1,1],[1,-2,-3]});T1 = T1/max(abs(T1(:)));
T2 = ncon({Tl,C},{[-1,-2,1],[1,-3]});T2 = T2/max(abs(T2(:)));

ncon({Tl,Tl},{[1,2,-1],[1,2,-2]});
delta = sqrt(ncon({T1-T2,conj(T1-T2)},{[1,2,3],[1,2,3]})); 

T3 = ncon({C,C},{[1,-1],[1,-2]});T3 = T3/max(abs(T3(:)));
T4 = ncon({C,C,T,T},{[1,2],[1,3],[2,4,-2],[3,4,-1]});T4 = T4/max(abs(T4(:)));

delta = sqrt(ncon({T3-T4,conj(T3-T4)},{[1,2,3],[1,2,3]}))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Compare tensors of the two different algorithms.
q = 2;
X = 4;
tol1 = 10^(-4);
tol2 = 10^(-4);

t = 1.2;

delta_4D = zeros(q,q,q,q);for i=1:q;delta_4D(i,i,i,i)=1;end
spin1 = zeros(q,q,q,q); spin1(1,1,1,1) = 1;
Qsq = sqrtm(ones(q)+(exp(1/t)-1)*eye(q));
A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
B = ncon({spin1,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
T = rand(X,q,X);T = T + permute(T,[3,2,1]);
C = rand(X,X);C = C + permute(C,[2,1]);
    
tic
[T1,C1,iter] = Potts_CTM2(T,C,B,A,q,X,tol1,t);
toc
tic
[T2,C2,iter2] = Potts_Alg(T,C,B,A,q,X,tol2,t);
toc

%T1./max(abs(T1(:)))
%T2./max(abs(T2(:)))
C1./max(abs(C1(:)))
C2./max(abs(C2(:)))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Plot magnetization vs temperature and magnetic field strength.
q = 2;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
spin1 = zeros(q,q,q,q);spin1(1,1,1,1)=1;

t_begin = 0.1;
t_step = 0.05;
t_end = 5;

h_begin = 0;
h_step = 0.05;
h_end = 4.9;

ts = t_begin:t_step:t_end;
hs = h_begin:h_step:h_end;

m = zeros(floor((t_end-t_begin)/t_step)+1,floor((h_end-h_begin)/h_step)+1);

X = [2];
tol = 10^(-4);

tic
for h = 1:numel(hs)
    for t = 1:numel(ts)
        
        Qsq = sqrtm([exp(1/ts(t)+hs(h)/ts(t)),1;1,exp(1/ts(t)-hs(h)/ts(t))]);
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        B = ncon({spin1,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
        [T,C,iter] = Potts_Alg(B,A,q,X,tol,t);
        T = permute(T,[2,1,3]);
        
        n = collapse(C,T,B)/collapse(C,T,A);
        m(h,t) = (q*n-1)/(q-1);
    end
end
toc

surf(ts,hs,m);
title('Magnetization as a function of temperature and magnetic field strenght for the 2D square Ising model')
xlabel('temperature')
ylabel('magnetic field strength')
zlabel('magnetization')
%}

function y = corrlen(T)
    %M = ncon({T,T},{[-1,1,-3],[-2,1,-4]});
    %M = reshape(M,size(T,1)^2,size(T,1)^2);
    %vals = eigs(M,2,'LM');
    %y = 1/abs(log(abs(vals(1))/abs(vals(2))));
    
    X = size(T,1);
    [~,vals] = eigs(@(x)mult1(T,T,x),X^2,2,'LM');
    y = 1/abs(log(abs(vals(1,1))/abs(vals(2,2))));
end

function y = mult1(M1,M2,C)
    si = size(M1);
    C = reshape(C,si(1),si(1));
    y = ncon({M1,M2,C},{[1,3,-2],[2,3,-1],[2,1]});
    y = reshape(y,si(1)^2,1);
end

function [C0,T0] = beginmatrices(Qsq,A,X)
    q = size(A,1);
    spin1_2D = zeros(q,q);spin1_2D(1,1)=1;
    spin1_3D = zeros(q,q,q);spin1_3D(1,1,1)=1;
    C0 = ncon({spin1_2D,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T0 = ncon({spin1_3D,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    
    while size(T0,1) < X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TA = ncon({T0,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TA},{[-1,1,2],[1,2,-2,-3,-4]});
        C0 = reshape(M,q*size(T0,1),q*size(T0,1));
        T0 = reshape(TA,[q*size(T0,1),q,q*size(T0,1)]);
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
    
    if size(T0,1) > X
        [U,~,~] = svd(C0);
        U_til = U(:,1:X);
        C0 = ncon({C0,U_til,U_til},{[1,2],[1,-1],[2,-2]});
        T0 = ncon({T0,U_til,U_til},{[1,-2,2],[1,-1],[2,-3]});
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
end