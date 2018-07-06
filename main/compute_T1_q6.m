set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
load('C_mxy_corr_q6X20-100tol6_FPCM_l.mat')

tcs = zeros(1,numel(X));
corrs = zeros(1,numel(X));

indices = [13,20,18,23,28,20,24,30,38,44,46,52,54,56,61,63,62,64,66,68,70,71,71,74,79,81,85,84,84,85,87,89,90,92,93,94,96,99,102,103, 102,104,104,107,107,108,109,110,113,115,113,115,116,117,120,124,125,128,123,124,127,126,127,127,130,128,128,130,130,130,131,132,132,133,134,135,136,137,136,137,139];
for x = 1:numel(X)
    index = indices(x);
    tcs(x) = ts(index);
    corrs(x) = corr2(x,index);
end

tcs = tcs(1:end);
corrs = corrs(1:end);

datax = corrs; datax = reshape(datax,numel(datax),1);
datay = tcs; datay = reshape(datay,numel(datay),1);
myfit = fit(datax,datay,fittype('a+b*log(c*x)^(-2)'),'StartPoint',[0.7,3,1]);
coeff = coeffvalues(myfit);Tc=coeff(1);a=coeff(2);b=coeff(3);
err_coeff = confint(myfit);err_Tc = err_coeff(2);

disp(['Tc from fit: ',num2str(Tc)]);
disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

subplot(2,1,1);
dataxshow = 0:0.001:log(b*datax(1))^(-2);
plot(log(b*datax).^(-2),datay,'bo');hold on;
plot(dataxshow,Tc+a*dataxshow,'r-','LineWidth',1);
set(gca,'fontsize',15)
legend('data','fit','Location','northeast');
xlabel('$1/l$','fontsize',25);
ylabel('$T^*_1$','fontsize',25);
