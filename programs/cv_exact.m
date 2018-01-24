function y = cv_exact(ts,delta_t)
    
    y = zeros(numel(ts),1);
    
    for t = 1:numel(ts)
        
        temp = ts(t)-delta_t;
        f0 = -(log(2)/2 + integral(@(x)integrand(x,temp),0,pi)/(2*pi))*temp;
        temp = ts(t);
        f1 = -(log(2)/2 + integral(@(x)integrand(x,temp),0,pi)/(2*pi))*temp;
        temp = ts(t)+delta_t;
        f2 = -(log(2)/2 + integral(@(x)integrand(x,temp),0,pi)/(2*pi))*temp;
        
        y(t) = -ts(t)*(f0-2*f1+f2)/(delta_t^2);
    
    end
end

function y = integrand(x,t)
    k = sinh(2/t)^(-2);
    y = log(((cosh(2/t))^2)+(sqrt(1+k^2-2*k*cos(2*x)))/k);
end