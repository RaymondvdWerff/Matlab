load('C_corr_q6X20-100tol6_FPCM_l.mat');

beginning = 1;
endings_ts = 20:39;
Tcs = zeros(1,numel(endings_ts));
errs = zeros(1,numel(endings_ts));
        
for ending_ts = 30:30   

    corr = zeros(1,numel(ts));
    cerr = zeros(1,numel(ts));

        for index = 1:numel(ts)
            xdata = xs(:,index);xdata = reshape(xdata,numel(xdata),1);
            ydata = ys(:,index);ydata = reshape(ydata,numel(ydata),1);
            myfit = fit(xdata,ydata,'poly1');
            coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);
            err_coeff = confint(myfit);corr_up = 1/err_coeff(3);corr_down = 1/err_coeff(4);

            corr(index) = (corr_up+corr_down)/2;
            cerr(index) = (corr_up-corr_down)/2;

        %         if corr(index) < 0
        %             ending_ts = index-1;
        %             break
        %         end

%             disp(['2nd ev from fit: ',num2str(b)]);
%             disp(['Error from fit: ',num2str(abs(b-err_b))]);
% 
%             xshow = 0:xdata(1)/100:xdata(1);
%             plot(xdata,ydata,'.');hold on;
%             plot(xshow,a*xshow+b);
        end
        
%         errorbar(ts(1:ending_ts),corr(1:ending_ts),cerr(1:ending_ts),'.');hold on;
        w = 1./cerr.^2;

        xdata = ts(1:ending_ts);xdata = reshape(xdata,numel(xdata),1);
        ydata = corr(1:ending_ts);ydata = reshape(ydata,numel(ydata),1);
        xdata_show = xdata(1):(xdata(2)-xdata(1))/100:xdata(end);

        myfittype = fittype('a*exp(b*(abs(c-x)/c)^(-1/2))');
        myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[0.01,4,0.9],'Weight',w(1:ending_ts));
        myfit = fit(xdata,ydata,myfittype,myfittopts);
        coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);Tc=coeff(3);
        err_coeff = confint(myfit);err_Tc = err_coeff(6);

        Tcs(ending_ts-19) = Tc;
        errs(ending_ts-19) = abs(Tc-err_Tc);
        
        disp(['Tc from fit: ',num2str(Tc)]);
        disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

        errorbar(xdata,ydata,cerr(1:ending_ts),'.');hold on;
        plot(xdata_show,a*exp(b*abs((xdata_show-Tc)/Tc).^(-1/2)));
        xlabel('temperature');
        ylabel('correlation length');
        title(['Correlation length as a function of temperature for the ',num2str(q),'-state clock model']);
        legend('data','fit');
end
%{
errorbar(endings_ts,Tcs,errs,'.');hold on;
xlabel('ending');
ylabel('Tc');
%}