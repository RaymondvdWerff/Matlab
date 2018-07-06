set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
load('C_mxy_corr_q6X20-100tol6_FPCM_r.mat')

tcs = zeros(1,numel(X));
corrs = zeros(1,numel(X));

indices = [347,313,304,292,275,267,256,256,243,221,205,190,179,180,176,170,167,166,165,163,164,158,154,141,130,124,123,119,114,113,108,109,108,107,101,91,89,85,82,77,72,70,69,66,61,61,61,58,56,54,56,54,53,53,53,53,53,50,52,48,47,45,44,39,35,33,28,27,27,26,26,23,21,22,20,18,16,13,12,13,12];
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
