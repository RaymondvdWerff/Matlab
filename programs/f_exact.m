function f = f_exact(t)
      
    f = -t*(log(2)/2 + integral(@(x)integrand(x,t),0,pi)/(2*pi))-1;
        
end

function y = integrand(x,t)
    k = sinh(1/t)^(-2);
    y = log(cosh(1/t)^2 + sqrt(1+k^2-2*k*cos(2*x))/k);
end