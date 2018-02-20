%Testing CTM and FPCM with the Potts/clock model.
%{
%Plot magnetization vs temperature.
q = 6;
X = [21,22,23];
tol1 = 6;
%tol2 = 4;
maxiter = 1000;

t_begin = 0.7;
t_step = 0.0001;
t_end = 0.98;

ts = t_begin:t_step:t_end;
%ts1 = 0.6:0.001:0.7; ts2 = 0.71:0.01:0.92; ts3 = 0.93:0.001:1.03; ts = [ts1,ts2,ts3];

for x = 1:numel(X)
    disp(['X = ' num2str(X(x))]);
    [mx1,my1,iters1,tictocs1] = compute_mxy(@Q_clock,q,X(x),10^(-tol1),maxiter,ts,@FPCM);
    save(['C_mxy_q' num2str(q) 'X' num2str(X(x)) 'tol' num2str(tol1) '_FPCM'],'ts','mx1','my1','iters1','tictocs1');
    %plot(ts,m1,'+-');hold on;
end

%X = 20:20:100;

%[m1,iters1,tictocs1] = compute_S(@Q_clock,q,X,10^(-tol1),maxiter,ts,@FPCM);
%[m2,iters2,tictocs2] = compute_f(@Q_clock,q,X,10^(-tol1),maxiter,ts,@FPCM);

%save(['C_Xm_q' num2str(q) 'X' num2str(X) 'tol' num2str(tol1) '_FPCM_bla'],'ts','m1','iters1','tictocs1');

%plot(ts,m1,'+-');hold on;
%plot(ts,my1,'+-');hold on;

%Exact solution Potts model:
%m = m_exact(ts);
%plot(ts,m,'+-');hold on;

%title(['Magnetization as a function of temperature for the ' num2str(q) '-state clock model (FPCM)']);
%xlabel('temperature');
%ylabel('magnetization');
%legend(['X = ' num2str(X(1))],['X = ' num2str(X(2))],['X = ' num2str(X(3))],['X = ' num2str(X(4))],['X = ' num2str(X(5))]);
%legend('X = 10','X = 20');

% subplot(2,2,1);
% plot(ts,mx1,'.');xlabel('T');ylabel('m_x');
% subplot(2,2,2);
% plot(ts,my1,'.');xlabel('T');ylabel('m_y');
% subplot(2,2,3);
% plot(ts,sqrt(mx1.^2+my1.^2),'.');xlabel('T');ylabel('|m|');
% subplot(2,2,4);
% plot3(mx1,my1,ts,'.');xlabel('m_x');ylabel('m_y');zlabel('T');

subplot(1,3,1);
plot(ts,mx1,'.');xlabel('T');ylabel('m_x');
axis([0.7 1 -0.8 0.8]);
subplot(1,3,2);
plot(ts,my1,'.');xlabel('T');ylabel('m_y');
subplot(1,3,3);
plot(ts,sqrt(mx1.^2+my1.^2),'.');xlabel('T');ylabel('|m|');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Plot iteratations vs bond dimension.
q = 6;
tol = 5;
maxiter = 1000;

temp = 0.8;
Qsq = sqrtm(Q_clock(q,temp,0));
delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end    
A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});

X = 10:2:100;
iters1 = zeros(1,numel(X));
iters2 = zeros(1,numel(X));

for x = 1:numel(X)
    disp(['X = ' num2str(X(x))]);
    [C,T] = beginmatrices(Qsq,A,X(x),1);
    [~,~,iter1] = CTM(A,C,T,X(x),10^(-tol),maxiter,temp);
    [~,~,iter2] = FPCM(A,C,T,X(x),10^(-tol),maxiter,temp);
    iters1(x) = iter1;
    iters2(x) = iter2;
end

plot(X,iters1,'+-');hold on;
plot(X,iters2,'x-');hold on;
title(['Iterations as a function of bond dimension for the ' num2str(q) '-state clock model at T = ' num2str(temp)]);
xlabel('bond dimension');
ylabel('number of iterations');
legend('CTM','FPCM');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine point of maximum slope of magnetization from fit
datax = ts; datax = datax(1:500);
datay = sqrt(mx1.^2+my1.^2); datay = datay(1:500);

P = polyfit(datax,datay,60);
D = polyder(P);

curve = polyval(D,datax);
[val,index] = min(curve);
ts(index)

% plot(datax,datay,'.-');hold on;
% plot(datax,polyval(P,datax),'.-');hold on;

plot(datax,polyval(D,datax),'.-');hold on;
plot(ts(index),curve(index),'o');
legend('fit','minimum');
%}
%{
X = 20:5:100;
tcs = zeros(1,numel(X));

for x = 1:numel(X)
    load(['C_Tc1_q6X', num2str(X(x)) ,'tol6_FPCM.mat']);
    datax = ts;
    datay = sqrt(mx1.^2+my1.^2);
    
    P = polyfit(datax,datay,6);
    D = polyder(P);

    % plot(datax,datay,'.-');hold on;
    % plot(datax,polyval(P,datax),'.-');hold on;

    curve = polyval(D,datax);
    [val,index] = min(curve);
    tcs(x) = ts(index);
end

%plot(1./log(X).^2,tcs,'o');hold on;
dataxfit = 1./log(X(6:end)).^2;
datayfit = tcs(6:end);
dataxfit = reshape(dataxfit,numel(dataxfit),1);
datayfit = reshape(datayfit,numel(datayfit),1);

myfit = fit(dataxfit,datayfit,fittype('a*x+b'));
coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
b
plot(myfit,dataxfit,datayfit);hold on;

% plot(datax(2:end),diff(datay)/(datax(2)-datax(1)),'.-');hold on;
% plot(datax,polyval(D,datax),'.-');hold on;
% plot(ts(index),curve(index),'o');
% legend('data','fit','minimum');
%}

%{
datax = ts(2:end);
datay = 1000*diff(sqrt(mx1.^2+my1.^2));

begin = 1;
ending = 100;
dataxfit = reshape(datax(begin:ending),ending-begin+1,1);
datayfit = reshape(datay(begin:ending),ending-begin+1,1);

dataxshow = datax(begin):0.0001:datax(ending);

myfittype = fittype('a*x^2+b*x+c');
myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1,1,1]);
myfit = fit(dataxfit,datayfit,'poly2');
%coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);c=coeff(3);

plot(datax,datay,'+-');hold on;
plot(myfit,dataxfit,datayfit);hold on;
%axis([ts(1) ts(end) 0 50])
legend('data','fit');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Compute Tc by divergence of correlation length
datax = X1;
datay = tcs1;

begin = 1;
ending = numel(X1);
dataxfit = reshape(datax(begin:ending),ending-begin+1,1);
datayfit = reshape(datay(begin:ending),ending-begin+1,1);

%dataxshow = datax(begin):0.001:datax(ending);
dataxshow = 20:0.1:90;

%myfittype = fittype('a*exp(b*(abs(xc-x)/xc)^(-1/2))');
myfittype = fittype('b/log(x/a)^2 + c');
myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[0.1,1,0.9]);
myfit = fit(dataxfit,datayfit,myfittype,myfittopts);
coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);Tc=coeff(3);
err_coeff = confint(myfit);err_Tc = err_coeff(5);

disp(['Tc from fit: ',num2str(Tc)]);
disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

plot(datax,datay,'o');hold on;
plot(dataxshow,b./log(dataxshow/a).^2 + Tc);
%plot(dataxshow,a*exp(b*(abs(Tc-dataxshow)/Tc).^(-1/2)),'x-');hold on;
%axis([ts(1) ts(end) 0 50])
legend('data','fit');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine magnetization as a function of iteration number/time
q = 6;
X = 20;
%tc = 1/log(1+sqrt(2));
%t = tc + 1e-3;
temp = 0.93;
maxiter1 = 500;
maxiter2 = 200;
tol = 1e-6;
tols = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];

for i = 1:10
    [m1,iters1,tictocs1,imarkers1,tmarkers1] = converge_m_FPCM(@Q_clock,q,X,tol,maxiter1,temp,tols);
    plot(tictocs1,m1,'.-');hold on;
end

%[m1,iters1,tictocs1,imarkers1,tmarkers1] = converge_m_CTM(@Q_clock,q,X,tol,maxiter1,temp,tols);
%[m2,iters2,tictocs2,imarkers2,tmarkers2] = converge_m_FPCM(@Q_clock,q,X,tol,maxiter2,temp,tols);
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
%{
%As a function of time.
plot(tictocs1,m1,'bo-');hold on;
%plot(tictocs1(end),m1(end),'*b');hold on;
plot(tictocs2,m2,'ro-');hold on;
%plot(tictocs2(end),m2(end),'*r');hold on;
for i = 1:numel(tmarkers1)
    plot([tmarkers1(i),tmarkers1(i)],[min(min(m1),min(m2)),max(max(m1),max(m2))],'b--');hold on;
end
for i = 1:numel(tmarkers2)
    plot([tmarkers2(i),tmarkers2(i)],[min(min(m1),min(m2)),max(max(m1),max(m2))],'r-.');hold on;
end
xlabel('time (s)');
ylabel('Magnetization');
title(['Magnetization as a function of time for the ' num2str(q) '-state Potts model for X = ' num2str(X) ' at T = ' num2str(temp)]);
legend('CTM','FPCM');
%}
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine the critical temperature for different bond dimensions.
format long

q = 2;
X = [10];
tol = 1e-3;
maxiter = 5000;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
tcs = zeros(numel(X),1);

for x = 1:numel(X)
    disp(['X = ' num2str(X(x))]);
    
    t_begin = 1.1345;
    t_step = 1e-4;
    t_end = 1.1349;

    tic
    for l = 1:2
        
        ts = t_begin:t_step:t_end;
        emptylist = zeros(numel(ts),1);
        m = emptylist;
        corr = emptylist;
        difference = 0;
        i = 1;
        counter = 0;
        
        for t = t_begin:t_step:t_end
            
            Qsq = sqrtm(Q_Potts(q,t));
            A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
            B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
            [C,T] = beginmatrices(Qsq,A,X(x));
            
            [C,T,~] = FPCM(A,C,T,q,X(x),tol,maxiter,t);

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
                    counter = 0;
                    tc = t;
                    index = i;
                    
                else
                    counter = counter + 1;
                    if counter > 2
                        break;
                    end
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

%Determine the critical temperature for different bond dimensions using interp.
%format long
%{
q = 2;
X = [20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100];
tol = 1e-3;
maxiter = 5000;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
tcs = zeros(numel(X),1);

for x = 1:numel(X)
    disp(['X = ' num2str(X(x))]);
    
    t_begin = 1.1345;
    t_step = 1e-4;
    t_end = 1.1349;
    
    i = 1;
    ts = t_begin:t_step:t_end;
    emptylist = zeros(numel(ts),1);
    
    corr = emptylist;
    
    for t = t_begin:t_step:t_end
        Qsq = sqrtm(ones(q)+(exp(1/t)-1)*eye(q));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        [C,T] = beginmatrices(Qsq,A,X(x));
    
        [C,T,~] = Potts_FPCM(A,C,T,q,X(x),tol,maxiter,t);
    
        corr(i) = corrlen(T);
        
        if i > 1
                if corr(i) > corr(i-1)
                    index = i;
                end
        end        
        i = i+1;
    end
    t_begin = ts(index)-t_step;
    t_end = ts(index)+t_step;
    t_step = t_step/10;
    
    ts = t_begin:t_step:t_end;
    [corr,~,~] = compute_corr(q,X(x),tol,maxiter,ts,@Potts_FPCM);
    
    t_step_int = t_step/1000;
    ts_int = t_begin:t_step_int:t_end;
    corr_int = interp1(ts,corr,ts_int,'spline');
    index = maximum(corr_int);
    tc = ts_int(index)
    tcs(x) = tc;
end
%}
%{
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
