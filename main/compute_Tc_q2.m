set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
format long
load('C_mxy_corr_q2X20-100tol8_CTM.mat')
load('C_mxy_corr_q2X20-100tol6_FPCM.mat')

tcs = zeros(1,numel(X));
corrs = zeros(1,numel(X));

for x = 1:numel(X)
    ts2 = ts(1):(ts(2)-ts(1))/10:ts(end);
    corr = interp1(ts,corr1(x,:),ts2,'spline');
    [val,index] = max(corr);
    tcs(x) = ts2(index);
    corrs(x) = val;    
end

tcs = tcs(1:end);
corrs = corrs(1:end);

datax = corrs; datax = reshape(datax,numel(datax),1);
datay = tcs; datay = reshape(datay,numel(datay),1);
dataxshow = 0:0.0001:0.002;
myfit = fit(datax,datay,fittype('a+b*x^(-c)'),'StartPoint',[2.269,0.34,1]);
coeff = coeffvalues(myfit);Tc=coeff(1);a=coeff(2);b=coeff(3);
err_coeff = confint(myfit);err_Tc = err_coeff(2);

disp(['Tc from fit: ',num2str(Tc)]);
disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

plot(1./datax,datay,'bo');hold on;
plot(dataxshow,Tc+a*dataxshow.^(b),'r-','LineWidth',1);
plot([0 0.002],[2/log(1+sqrt(2)),2/log(1+sqrt(2))],'k--');
set(gca,'fontsize',15)
legend('data','fit','exact T_c','Location','northwest');
xlabel('$1/\xi$','fontsize',25);
ylabel('$T_C^*$','fontsize',25);
title('\textbf{FPCM}');
