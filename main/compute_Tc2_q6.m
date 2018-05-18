load('C_corr_q6X20-100tol6_FPCM_r.mat');

q = 6;
beginning_X = 1;
beginnings_ts = 45:70;
Tcs = zeros(1,numel(beginnings_ts));
errs = zeros(1,numel(beginnings_ts));
        
for beginning_ts = 55:55  

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
        %             ending = index-1;
        %             break
        %         end

%             disp(['2nd ev from fit: ',num2str(b)]);
%             disp(['Error from fit: ',num2str(abs(b-err_coeff(4)))]);
% 
%             xshow = 0:xdata(1)/100:xdata(1);
%             plot(xdata,ydata,'.');hold on;
%             plot(xshow,a*xshow+b);
        end
        
%         errorbar(ts(beginning_ts:end),corr(beginning_ts:end),cerr(beginning_ts:end),'.');hold on;
        w = 1./cerr.^2;
       
        xdata = ts(beginning_ts:end);xdata = reshape(xdata,numel(xdata),1);
        ydata = corr(beginning_ts:end);ydata = reshape(ydata,numel(ydata),1);
        xdata_show = xdata(1):(xdata(2)-xdata(1))/100:xdata(end);
        
        myfittype = fittype('a*exp(b*(abs(c-x)/c)^(-1/2))');
        myfittopts = fitoptions('Method','NonLinearLeastSquares','StartPoint',[2,1,0.9],'Weight',w(beginning_ts:end));
        myfit = fit(xdata,ydata,myfittype,myfittopts);
        coeff = coeffvalues(myfit);a=coeff(1);b=coeff(2);Tc=coeff(3);
        err_coeff = confint(myfit);err_Tc = err_coeff(6);

%         Tcs(beginning_ts-beginnings_ts(1)+1) = Tc;
%         errs(beginning_ts-beginnings_ts(1)+1) = abs(Tc-err_Tc);
        
        disp(['Tc from fit: ',num2str(Tc)]);
        disp(['Error from fit: ',num2str(abs(Tc-err_Tc))]);

        errorbar(xdata,ydata,cerr(beginning_ts:end),'.');hold on;
        plot(xdata_show,a*exp(b*abs((xdata_show-Tc)/Tc).^(-1/2)));
        xlabel('temperature');
        ylabel('correlation length');
        title(['Correlation length as a function of temperature for the ',num2str(q),'-state clock model']);
        legend('data','fit');
end
%{
errorbar(beginnings_ts,Tcs,errs,'.');hold on;
xlabel('beginning');
ylabel('Tc');
%}