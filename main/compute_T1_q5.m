set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
load('C_mxy_corr_q5X20-100tol6_FPCM_l.mat')

tcs = zeros(1,numel(X));
corrs = zeros(1,numel(X));

indices = [20,26,26,33,31,34,27,30,35,41,47,48,47,50,53,56,58,57,63,64,65,68,68,68,67,71,74,73,73,74,75,78,79,80,81,83,84,86,87,87,87,88,87,88,88,89,92,94,94,96,99,99,100,101,101,102,103,103,103,104,106,107,107,108,109,109,108,109,110,110,110,111,112,112,113,113,114,115,115,115,116];      
for x = 1:numel(X)
    index = indices(x);
    tcs(x) = ts(index);
    corrs(x) = corr2(x,index);
end

tcs = tcs(1:end);
corrs = corrs(1:end);

datax = corrs; datax = reshape(datax,numel(datax),1);
datay = tcs; datay = reshape(datay,numel(datay),1);
myfit = fit(datax,datay,fittype('a+b*log(c*x)^(-2)'),'StartPoint',[0.9,3,1]);
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
