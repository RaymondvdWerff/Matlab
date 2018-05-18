%Testing CTM and FPCM with the Potts/clock model.

%{
%Plot OBSERVABLE vs temperature.
q = 6;
X = 3;
tol1 = 1e-6;
tol2 = 1e-3;
maxiter1 = 5000;
maxiter2 = 2000;

t_begin = 0.5;
t_step = 0.001;
t_end = 1.6;

ts = t_begin:t_step:t_end;

%[mx1,my1,iters1,tictocs1] = compute_m(@Q_clock,q,X,tol1,maxiter1,ts,@FPCM,0);
[m2,corr2,iters2,tictocs2] = compute_m_and_corr(@Q_clock,q,X,tol2,maxiter2,ts,@FPCM,0); 

plot(ts,corr2,'-');
%legend('CTM','FPCM');
%suptitle(['Magnetization of the ',num2str(q),'-state Clock model for \chi = ',num2str(X)]);

% mx = zeros(numel(X),numel(tol2),numel(ts));
% my = zeros(numel(X),numel(tol2),numel(ts));
% iters = zeros(numel(X),numel(ts));
% tictocs = zeros(numel(X),numel(ts));
% 
% for x = 1:numel(X)
%     disp(['X = ',num2str(X(x))]);
%     [mx1,my1,iters1,tictocs1] = compute_mxy(@Q_clock,q,X(x),tol2,maxiter2,ts,@FPCM,0);
%     mx(x,:,:) = mx1;
%     my(x,:,:) = my1;
%     iters(x,:) = iters1;
%     tictocs(x,:) = tictocs1;
% end

% for x = 1:numel(X)
%     mx1 = reshape(mx(x,1,:),1,numel(mx(x,1,:)));
%     my1 = reshape(my(x,1,:),1,numel(my(x,1,:)));
%     plot(mx1./sqrt(mx1.^2+my1.^2),my1./sqrt(mx1.^2+my1.^2),'.');hold on;
% end

% xlabel('m_x');
% ylabel('m_y');
% title(['Magnetization of the ',num2str(q),'-state clock model at T = ',num2str(temp)]);

%angles = 360*atan2(mx2,my2)/(2*pi);angles = angles(:);histogram(angles,360);hold on;
%xlim([-180,180]);xticks(-180:60:180);

% for i = 1:numel(tol2)
%     subplot(2,3,i);
%     plot(mx1(i,:),my1(i,:),'b.');hold on;
%     bnd = max(mx1(1,:));
%     plot([-bnd,bnd],[-bnd*tan(pi/6),bnd*tan(pi/6)],'k-.');
%     plot([bnd,-bnd],[-bnd*tan(pi/6),bnd*tan(pi/6)],'k-.');
%     plot([0,0],[-bnd,bnd],'k-.');
%     plot([-bnd/tan(pi/3),bnd/tan(pi/3)],[-bnd,bnd],'k:');
%     plot([-bnd/tan(pi/3),bnd/tan(pi/3)],[bnd,-bnd],'k:');
%     plot([-bnd,bnd],[0,0],'k:');
% %     plot(mx1(i,37),my1(i,37),'ro');
% %     plot(mx1(i,57),my1(i,57),'ro');
% %     plot(mx1(i,129),my1(i,129),'ro');
% %     plot(mx1(i,203),my1(i,203),'ro');
% %     plot(mx1(i,262),my1(i,262),'ro');
% %     plot(mx1(i,412),my1(i,412),'ro');
% %     plot(mx1(i,446),my1(i,446),'ro');
%     xlabel('m_x');
%     ylabel('m_y');
%     title(['tol = ', num2str(tol2(i))]);
% end

% for i = 1:numel(tol2)
%    subplot(2,3,i);
%    angles = 360*atan2(mx1(i,:),my1(i,:))/(2*pi);
%    histogram(angles,360);
%    xlim([-180,180]);
%    xticks(-180:60:180);
% end

%suptitle('Magnetization of the 6-state clock model at T = 0.8');

%save(['C_m_q' num2str(q) 'X' num2str(X) 'tol' num2str(tol1) '_FPCM'],'ts','mx1','my1','mx2','my2','iters2','tictocs2');

% subplot(2,2,1);
% plot(ts,mx1,'.');xlabel('T');ylabel('m_x');
% subplot(2,2,2);
% plot(ts,my1,'.');xlabel('T');ylabel('m_y');
% subplot(2,2,3);
% plot(ts,sqrt(mx1.^2+my1.^2),'.');xlabel('T');ylabel('|m|');
% subplot(2,2,4);
% plot3(mx1,my1,ts,'.');xlabel('m_x');ylabel('m_y');zlabel('T');

% subplot(1,3,1);
% plot(ts,mx2,'.-');xlabel('T');ylabel('m_x');
% subplot(1,3,2);
% plot(ts,my2,'.-');xlabel('T');ylabel('m_y');
% subplot(1,3,3);
% plot(ts,sqrt(mx2.^2+my2.^2),'.-');xlabel('T');ylabel('|m|');
% suptitle(['Magnetization of the ',num2str(q),'-state Clock model for \chi = ',num2str(X)]);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine magnetization as a function of iteration number/time
q = 2;
X = 600;
%tc = 1/log(1+sqrt(2));
%t = tc + 1e-3;
%temp = (1+0.001)*2/log(1+sqrt(2));
%temp = (1+0.0001)*2/log(1+sqrt(2));

%temp = 0.825;
for i = 1:numel(temp)
    disp(['i = ',num2str(i)]);
    maxiter1 = 30000;
    maxiter2 = 600;
    tol1 = 1e-6;
    tol2 = 1e-6;
    tols = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];

    % for i = 10:10:100
    %     [m1,iters1,tictocs1,imarkers1,tmarkers1] = converge_m_FPCM(@Q_clock,q,X,tol,maxiter1,temp,tols,i);
    %     plot(tictocs1(end),m1(end),'o');hold on;
    % end

    %[m1,sv1,iters1,tictocs1,imarkers1,tmarkers1] = converge_m_CTM(@Q_clock,q,X,tol1,maxiter1,temp(i),tols);
    %[m2,sv2,iters2,tictocs2,imarkers2,tmarkers2] = converge_m_FPCM(@Q_clock,q,X,tol2,maxiter2,temp(i),tols);
    %[m3,sv3,iters3,tictocs3,imarkers3,tmarkers3] = converge_m_CTM(@Q_clock,q,X,tol,maxiter1,temp2,tols);
    %[m4,sv4,iters4,tictocs4,imarkers4,tmarkers4] = converge_m_FPCM(@Q_clock,q,X,tol,maxiter2,temp2,tols);

    % subplot(1,2,1);
    % semilogy(tictocs1,m1,'.-');hold on;
    % semilogy(tictocs2,m2,'.-');
    % xlabel('time (s)');ylabel('magnetization');legend('CTM','FPCM');
    % title('T/Tc-1 = 0.001');
    % 
    % subplot(1,2,2);
    % semilogy(tictocs3,m3,'.-');hold on;
    % semilogy(tictocs4,m4,'.-');
    % xlabel('time (s)');ylabel('magnetization');legend('CTM','FPCM');
    % title('T/Tc-1 = 0.0001');

%     subplot(2,1,1);
%     semilogy(tictocs1,m1,'.-');hold on;
%     semilogy(tictocs2,m2,'.-');
%     xlabel('time (s)');ylabel('magnetization');legend('CTM','FPCM');
%     
    % subplot(1,3,2);
    % semilogy(tictocs1,abs(m1-m1(end)),'.-');hold on;
    % semilogy(tictocs2,abs(m2-m2(end)),'.-');
    % xlabel('time (s)');ylabel('difference of magnetization');legend('CTM','FPCM');
    % 
    %subplot(2,3,i);
    semilogy(tictocs1/3600,m1,'.-');hold on;
    semilogy(tictocs2/3600,m2,'.-');hold on;
    title(['temp = (1+0.001)*Tc = ',num2str(temp(i))]);
    xlabel('time (h)');ylabel('error of magnetization');legend('CTM','FPCM');
end
%suptitle(['Convergence of the magnetization and singular value spectrum for the ' num2str(q) '-state Clock model for \chi = ' num2str(X)]);

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
%Determination of critical temperature from fitting data
q = 6;
tol1 = 6;

load('C_Tc2_q6X20-100tol6_FPCM_b');
load(['C_corr2_q' num2str(q) 'X' num2str(X(1)) 'tol' num2str(tol1) '_FPCM']);

X = X(1:end);
tcs = tcs(1:end);
corrs = zeros(1,numel(X));
index = 1;
%index = 65;
ts(index)

for x = 1:numel(X)
    load(['C_corr2_q' num2str(q) 'X' num2str(X(x)) 'tol' num2str(tol1) '_FPCM']);
    corrs(x) = corr1(index);
end

dataxfit = reshape(corrs,numel(corrs),1);
datayfit = reshape(tcs,numel(tcs),1); 

myfit = fit(dataxfit,datayfit,fittype('Tc+a/log(b*x)^2'),'StartPoint',[0.7,-20,1000]);
coeff = coeffvalues(myfit);Tc=coeff(1);a=coeff(2);b=coeff(3);
err_coeff = confint(myfit);err_Tc = err_coeff(2);

disp(['Tc from fit: ',num2str(Tc)]);
disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

dataxshow = dataxfit(1):0.01:dataxfit(end);

plot(1./log(b*dataxfit).^2,datayfit,'o');hold on;
plot(1./log(b*dataxshow).^2,Tc+a./log(b*dataxshow).^2);hold on;
legend('data','fit');
%}

%{
for begin = 1:3
    load('C_Tc1_q6X20-100tol6_FPCM_b');
    load('C_corr_q6X20tol6_FPCM');
    X = X(begin:end);
    tcs = tcs(begin:end);

    tcin = zeros(1,numel(ts));
    tcout = zeros(1,numel(ts));

    for index = 1:numel(ts)

        corrs = zeros(1,numel(X));
        for x = 1:numel(X)
           load(['C_corr_q' num2str(q) 'X' num2str(X(x)) 'tol' num2str(tol1) '_FPCM']);
           corrs(x) = corr1(index);
        end

        dataxfit = reshape(corrs,numel(corrs),1);
        datayfit = reshape(tcs,numel(tcs),1); 

        myfit = fit(dataxfit,datayfit,fittype('Tc+a/log(b*x)^2'),'StartPoint',[0.8,-20,1000]);
        coeff = coeffvalues(myfit);Tc=coeff(1);a=coeff(2);b=coeff(3);
        err_coeff = confint(myfit);err_Tc = err_coeff(2);

        tcin(index) = ts(index);
        tcout(index) = Tc;
    end

    plot(tcin,tcout,'.');hold on;
end
plot([ts(1),ts(end)],[ts(1),ts(end)]);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Plot iteratations vs bond dimension.
q = 6;
tol = 6;
maxiter = 1200;

temp = 0.8;
Qsq = sqrtm(Q_clock(q,temp,0));
delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end    
A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
Qsq2 = sqrtm(Q_clock(q,temp,1e-4));
A2 = ncon({delta_4D,Qsq2,Qsq2,Qsq2,Qsq2},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
       
X = 10:1:100;
%iters1 = zeros(1,numel(X));
iters2 = zeros(1,numel(X));

for x = 1:numel(X)
    disp(['X = ' num2str(X(x))]);
    [C,T] = beginmatrices(Qsq2,A2,X(x),0);
    %[~,~,iter1] = CTM(A,C,T,X(x),10^(-tol),maxiter,temp);
    [~,~,iter2] = FPCM(A,C,T,X(x),10^(-tol),maxiter,temp);
    %iters1(x) = iter1;
    iters2(x) = iter2;
end

%plot(X,iters1,'+-');hold on;
plot(X,iters2,'x-');hold on;
title(['Iterations as a function of bond dimension for the ' num2str(q) '-state clock model at T = ' num2str(temp)]);
xlabel('bond dimension');
ylabel('number of iterations');
%legend('CTM','FPCM');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Determine point of maximum slope of magnetization from fit
beginning = 1;
ending = numel(ts);
datax = ts;datax = datax(beginning:1:ending);
datay = sqrt(mx1.^2+my1.^2); datay = datay(beginning:1:ending);
dataxshow = datax(1):0.00001:datax(end);

P = polyfit(datax,datay,4);
D = polyder(P);

curve1 = polyval(P,dataxshow);
curve2 = polyval(D,dataxshow);
[val,index] = min(curve2);
dataxshow(index)

subplot(1,2,1);
plot(datax,datay,'.-');hold on;
plot(dataxshow,curve1,'-');hold on;
plot(dataxshow(index),curve1(index),'o');
legend('data','fit','steepest slope');

subplot(1,2,2);
plot(datax(2:end)-(datax(2)-datax(1))/2,diff(datay)/(datax(2)-datax(1)),'.-');hold on;
plot(dataxshow,curve2,'-');hold on;
plot(dataxshow(index),curve2(index),'o');
legend('data','fit','minimum');
%}
%{
X = 20:5:100;
tcs = zeros(1,numel(X));

for x = 1:numel(X)
    load(['C_Tc1_q6X', num2str(X(x)) ,'tol6_FPCM.mat']);
    ending = numel(ts)-25;
    ts = ts(1:ending); mx1 = mx1(1:ending); my1 = my1(1:ending);
    y = sqrt(mx1.^2+my1.^2);
    datax = ts(1:10:end);
    datay = y(1:10:end);
    
    P = polyfit(datax,datay,4);
    D = polyder(P);
    
    %plot(datax,datay,'.-');hold on;
    plot(ts(2:end),10000*diff(y),'.-');hold on;
    plot(ts,polyval(D,ts),'.-');hold on;

    curve = polyval(D,ts);
    [~,index] = min(curve);
    tcs(x) = ts(index);
end

%plot(1./log(X).^2,tcs,'o');hold on;
% dataxfit = 1./log(X(1:end)).^2;
% datayfit = tcs(1:end);
% dataxfit = reshape(dataxfit,numel(dataxfit),1);
% datayfit = reshape(datayfit,numel(datayfit),1);
% 
% myfit = fit(dataxfit,datayfit,fittype('a*x+b'));
% coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
% b
% plot(myfit,dataxfit,datayfit);hold on;

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
%Determination of critical temperature by extrapolating correlation length
X_beg = 1;
ts_beg = 167;
ts_end = 231;
ts_fit = ts_beg:ts_end;

corr = zeros(1,numel(ts(ts_beg:ts_end)));
for index = ts_fit
    xdata = xs(X_beg:end,index);xdata = reshape(xdata,numel(xdata),1);
    ydata = ys(X_beg:end,index);ydata = reshape(ydata,numel(ydata),1);
    myfit = fit(xdata,ydata,'poly1');
    coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
    err_coeff = confint(myfit);err_b = err_coeff(4);

    corr(index-ts_beg+1) = 1/b;
%     disp(['2nd ev from fit: ',num2str(b)]);
%     disp(['Error from fit: ',num2str(abs(b-err_b))]);
% 
%     xshow = 0:xdata(1)/100:xdata(1);
%     plot(xdata,ydata,'.');hold on;
%     plot(xshow,a*xshow+b);
end

% for i = 1:numel(X)
%     plot(ts,1./ys(i,:),'.-');hold on
% end
%plot(ts,corr,'.-');hold on;

Tc = 0.91:0.001:0.965;
Tc = 0.958:0.958;
errs_a = zeros(1,numel(Tc));
errs_b = zeros(1,numel(Tc));

for x = 1:numel(Tc)
    xdata = (abs(ts(ts_beg:ts_end)-Tc(x))/Tc(x)).^(-1/2);xdata = reshape(xdata,numel(xdata),1);
    ydata = log(corr);ydata = reshape(ydata,numel(ydata),1);
    
    myfittype = fittype('a+x*b');
    myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1,1]);
    myfit = fit(xdata,ydata,myfittype,myfittopts);
    coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
    err_coeff = confint(myfit);err_a = err_coeff(2);err_b = err_coeff(4);
    errs_a(x) = abs(err_a-a);
    errs_b(x) = abs(err_b-b);
    
%     plot(xdata,ydata,'.');hold on;
%     plot(xdata,a+xdata*b);
    
end

plot(ts(ts_beg:ts_end),corr,'.');hold on;
plot(ts(ts_beg:ts_end),exp(a)*exp(b*(abs(ts(ts_beg:ts_end)-Tc)/Tc).^(-1/2)))

% subplot(1,2,1);
% plot(Tc,errs_a,'.');hold on;
% subplot(1,2,2);
% plot(Tc,errs_b,'.');hold on;
%  
% [~,index1] = min(errs_a);
% [~,index2] = min(errs_b);
% Tc = Tc(index1);
% disp(['Tc = ',num2str(Tc)]);
%}
%{
xdata = (abs(ts(beginning_ts:end)-Tc)/Tc).^(-1/2);xdata = reshape(xdata,numel(xdata),1);
ydata = log(corr(beginning_ts:end));ydata = reshape(ydata,numel(ydata),1);

myfittype = fittype('a+x*b');
myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1,1]);
myfit = fit(xdata,ydata,myfittype,myfittopts);
coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);

xdata = ts(beginning_ts:end);xdata = reshape(xdata,numel(xdata),1);
ydata = corr(beginning_ts:end);ydata = reshape(ydata,numel(ydata),1);
myfittype = fittype('a*exp(b*(abs(c-x)/c)^(-1/2))');
myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[exp(a),b,Tc]);
myfit = fit(xdata,ydata,myfittype,myfittopts);
coeff = coeffvalues(myfit);a2=coeff(1);b2=coeff(2);Tc2=coeff(3);
err_coeff = confint(myfit);err_Tc2 = err_coeff(6);

disp(['Tc from fit: ',num2str(Tc2)]);
disp(['Error from fit: ',num2str(abs(Tc2-err_Tc2))]);

plot(xdata,ydata,'.');hold on;
%plot(xdata,exp(a)*exp(b*abs((xdata-Tc)/Tc).^(-1/2)));
plot(xdata,a2*exp(b2*abs((xdata-Tc2)/Tc2).^(-1/2)));
%axis([min(ts) max(ts) min(corr) max(corr)])
legend('data','fit');

%}
%{
beginnings = 17:17;
tcs = zeros(1,numel(beginnings));

for j = 1:numel(beginnings)
    corr = zeros(1,numel(ts));
    for index = 1:numel(ts)
        xdata = xs(beginnings(j):end,index);xdata = reshape(xdata,numel(xdata),1);
        ydata = ys(beginnings(j):end,index);ydata = reshape(ydata,numel(ydata),1);
        myfit = fit(xdata,ydata,'poly1');
        coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
        err_coeff = confint(myfit);err_a = err_coeff(2);err_b = err_coeff(4);

        corr(index) = 1/b;
%         disp(['2nd ev from fit: ',num2str(b)]);
%         disp(['Error from fit: ',num2str(abs(b-err_b))]);
% 
%         xshow = 0:xdata(1)/100:xdata(1);
%         plot(xdata,ydata,'.');hold on;
%         plot(xshow,a*xshow+b);
    end
    %plot(ts,corr,'.');hold on;
    
    endings = 20:20;
    Tcs = zeros(1,numel(endings));
    for i = 1:numel(endings)
        ending = endings(i);
        % xdata = ts(1:ending);xdata = reshape(xdata,numel(xdata),1);
        % ydata = corr(1:ending);ydata = reshape(ydata,numel(ydata),1);

        % myfittype = fittype('a*exp(b*(abs(Tc-x)/Tc)^(-1/2))');
        % myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[0.04,2.05,0.706]);
        % myfit = fit(xdata,ydata,myfittype,myfittopts);
        % coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);Tc=coeff(3);
        % err_coeff = confint(myfit);err_Tc = err_coeff(5);
        % 
        % plot(xdata,ydata,'.');hold on;
        % plot(xdata,c*exp(d*(abs(Td-xdata)/Td).^(-1/2)));
        % plot(xdata,a*exp(b*(abs(Tc-xdata)/Tc).^(-1/2)));
        %axis([xdata(1) xdata(end) 0 ydata(end)]);
        %plot(xdata,0.0028*exp(3.055*(abs(0.71827-xdata)/0.71827).^(-1/2)));

        Tc = 0.68:0.0001:0.76;
        errs_a = zeros(1,numel(Tc));
        errs_b = zeros(1,numel(Tc));

        for x = 1:numel(Tc)
            xdata = (abs(ts(1:ending)-Tc(x))/Tc(x)).^(-1/2);xdata = reshape(xdata,numel(xdata),1);
            ydata = log(corr(1:ending));ydata = reshape(ydata,numel(ydata),1);

            myfittype = fittype('a+x*b');
            myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1,1]);
            myfit = fit(xdata,ydata,myfittype,myfittopts);
            coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
            err_coeff = confint(myfit);err_a = err_coeff(2);err_b = err_coeff(4);
            errs_a(x) = abs(err_a-a);
            errs_b(x) = abs(err_b-b);
        end

    %     subplot(1,2,1);
    %     plot(Tc,errs_a,'.');hold on;
    %     plot(Tc,errs_b,'.');hold on;

        [~,index1] = min(errs_a);
        [~,index2] = min(errs_b); 
        %Tc(index1)
        Tcs(i) = Tc(index1);
    end

    tcs(j) = min(Tcs);
end

plot(endings,Tcs,'.');hold on;

% xdata = (abs(ts(beginning:ending)-Tc)/Tc).^(-1/2);xdata = reshape(xdata,numel(xdata),1);
% ydata = log(corr(beginning:ending));ydata = reshape(ydata,numel(ydata),1);
% 
% myfittype = fittype('a+x*b');
% myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[1,1]);
% myfit = fit(xdata,ydata,myfittype,myfittopts);
% coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
% err_coeff = confint(myfit);err_a = err_coeff(2);
% errs_a(x) = abs(err_a-a);
% 
% subplot(1,2,2);
% plot(xdata,ydata,'.');hold on;
% plot(xdata,a+xdata*b);hold on;

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
