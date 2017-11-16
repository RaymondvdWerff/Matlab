N = 10000
t = [0:N]/N;
y1 = sin(2*pi*t);
y2 = exp(-5*t);
y3 = cos(t).*exp(-t);

integral_1 = sum((y1(1:N)+y1(2:N+1)).*(diff(t)))/2
integral_2 = sum((y2(1:N)+y2(2:N+1)).*(diff(t)))/2
integral_3 = sum((y3(1:N)+y3(2:N+1)).*(diff(t)))/2
