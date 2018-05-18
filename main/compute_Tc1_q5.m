load('C_corr_q5X20-100tol6_FPCM.mat');

i=1;
% ts_ends = 20:1:100;
% Tcs = zeros(1,numel(ts_ends));
% errs = zeros(1,numel(ts_ends));
       
for bla = 1:1

    ts_beg = 67;
    ts_end = 101;
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
    xdata_show = ts(ts_beg):(ts(2)-ts(1))/100:ts(ts_end+3);

    myfittype = fittype('a*exp(b*(abs(c-x)/c)^(-1/2))');
    myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[0.01,2,1],'Weight',w);
    myfit = fit(xdata,ydata,myfittype,myfittopts);
    coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);Tc=coeff(3);
    err_coeff = confint(myfit);err_Tc = err_coeff(6);

%     Tcs(i) = Tc;
%     errs(i) = abs(Tc-err_Tc);
% 
    disp(['Tc from fit: ',num2str(Tc)]);
    disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

    errorbar(xdata,ydata,cerr,'b.');hold on;
    plot(xdata_show,a*exp(b*abs((xdata_show-Tc)/Tc).^(-1/2)),'r-');hold on;
    plot([Tc,Tc],[0,3000],'k--');hold on;
    xlabel('temperature');
    ylabel('correlation length');
    title('Correlation length as a function of temperature for the 5-state clock model');
    legend('data','fit');
    i = i+1;
end

% errorbar(ts_ends,Tcs,errs,'.');hold on;
% xlabel('ending');
% ylabel('Tc');
