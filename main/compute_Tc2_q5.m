load('C_evs_q5X20-100tol6_FPCM_r.mat');

xs = log(evs(:,:,2)./evs(:,:,5));
ys = -log(evs(:,:,2)); 
        
i=1;
ts_ends = ts(41:101);
Tcs = zeros(1,numel(ts_ends));
errs = zeros(1,numel(ts_ends));
       
for ts_end = 41:101
    
    ts_beg = 21;
%     ts_end = 101;
    ts_fit = ts_beg:ts_end;
    
    corr = zeros(1,numel(ts_fit));
    cerr = zeros(1,numel(ts_fit));
    
    for index = ts_fit
        xdata = xs(:,index);xdata = reshape(xdata,numel(xdata),1);
        ydata = ys(:,index);ydata = reshape(ydata,numel(ydata),1);
        myfit = fit(xdata,ydata,'poly1');
        err_coeff = confint(myfit);corr_up = 1/err_coeff(3);corr_down = 1/err_coeff(4);

        corr(index-ts_beg+1) = (corr_up+corr_down)/2;
        cerr(index-ts_beg+1) = (corr_up-corr_down)/2;

%         coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);err_b=err_coeff(3);
%         disp(['2nd ev from fit: ',num2str(b)]);
%         disp(['Error from fit: ',num2str(abs(b-err_b))]);
% 
%         xshow = 0:xdata(1)/100:xdata(1);
%         plot(xdata,ydata,'.');hold on;
%         plot(xshow,a*xshow+b);
    end
        
%     plot(ts(ts_beg:ts_end),corr,'.-');hold on;
    w = 1./cerr.^2;

    xdata = ts(ts_beg:ts_end);xdata = reshape(xdata,numel(xdata),1);
    ydata = corr;ydata = reshape(ydata,numel(ydata),1);
    xdata_show = ts(ts_beg)-3*(ts(2)-ts(1)):(ts(2)-ts(1))/100:ts(ts_end);

    myfittype = fittype('a*exp(b*(abs(c-x)/c)^(-1/2))');
    myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[0.01,4,0.9],'Weight',w);
    myfit = fit(xdata,ydata,myfittype,myfittopts);
    coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);Tc=coeff(3);
    err_coeff = confint(myfit);err_Tc = err_coeff(6);

    Tcs(i) = Tc;
    errs(i) = abs(Tc-err_Tc);

%     disp(['Tc from fit: ',num2str(Tc)]);
%     disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

%     errorbar(xdata,ydata,cerr,'b.');hold on;
%     plot(xdata_show,a*exp(b*abs((xdata_show-Tc)/Tc).^(-1/2)),'r-');hold on;
%     plot([Tc,Tc],[0,max(a*exp(b*abs((xdata_show-Tc)/Tc).^(-1/2)))],'k--');hold on;
%     set(gca,'fontsize',15)
%     xlabel('$T$','fontsize',25);
%     ylabel('$\xi$','fontsize',25);
%     legend('data','fit');
    i = i+1;
end
% for x = 1:numel(X)
%     plot(ts,1./ys(x,:));hold on;
% end

subplot(2,1,2);
errorbar(ts_ends,Tcs,errs,'b.');hold on;
plot([ts_ends(1),ts_ends(end)],[max(Tcs+errs),max(Tcs+errs)],'k--');
plot([ts_ends(1),ts_ends(end)],[min(Tcs-errs),min(Tcs-errs)],'k--');
set(gca,'fontsize',15);
xlabel('$T_{end}$');
ylabel('$T_2$');
xlim([ts_ends(1),ts_ends(end)]);