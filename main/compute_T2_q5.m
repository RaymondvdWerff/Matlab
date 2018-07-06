set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
load('C_mxy_corr_q5X20-100tol6_FPCM_r.mat')

tcs = zeros(1,numel(X));
corrs = zeros(1,numel(X));

indices = [257,245,232,216,206,191,184,176,175,162,153,149,144,141,140,137,131,131,118,119,117,112,107,105,103,99,96,94,90,85,83,80,76,70,70,71,66,64,58,55,50,49,48,47,47,44,43,41,40,40,40,39,37,35,34,35,35,36,32,32,31,30,28,26,23,23,22,20,20,18,18,16,15,14,13,13,13,12,11,9,6];
for x = 1:numel(X)
%     y = diff(sqrt(mx2(x,:).^2+my2(x,:).^2));
%     [~,index] = min(y);
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

subplot(2,1,2);
dataxshow = 0:0.001:log(b*datax(1))^(-2);
plot(log(b*datax).^(-2),datay,'bo');hold on;
plot(dataxshow,Tc+a*dataxshow,'r-','LineWidth',1);
set(gca,'fontsize',15)
legend('data','fit','Location','northwest');
xlabel('$1/l$','fontsize',25);
ylabel('$T^*_2$','fontsize',25);
