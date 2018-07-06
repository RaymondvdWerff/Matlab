function f = f_exact(ts)

    f = zeros(1,numel(ts)); 
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        f(t) = -temp*(log(2)/2 + integral(@(x)integrand(x,temp),0,pi)/(2*pi));
    
    end    
end

function y = integrand(x,t)
    k = sinh(2/t)^(-2);
    y = log(cosh(2/t)^2 + sqrt(1+k^2-2*k*cos(2*x))/k);
end