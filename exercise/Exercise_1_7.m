for N = 10:1000
    t = [0:N]/N;
    y= exp(-t);
    number(N-9) = N;
    integral = sum((y(1:N)+y(2:N+1)).*(diff(t)))/2;
    exact = 1-exp(-1);
    error(N-9) = abs(integral-exact);
end

loglog(number,error)
