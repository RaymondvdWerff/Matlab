x = -1:0.1:2;
xint = -1:0.01:2;
y = sin(x);
y2 = sin(xint);

yint = interp1(x,y,x2,'spline');
plot(x,y);hold on;
plot(xint,y2,'o-');
plot(xint,yint,'x-');