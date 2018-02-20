
x1 = abs(mx1);
y1 = abs(my1);
x2 = x1;
y2 = y1;

x3 = abs(mx2);
y3 = abs(my2);
x4 = x3;
y4 = y3;

for i = 1:numel(ts)
    
    if atan(y1(i)/x1(i)) > pi/3
        x2(i) = 0.5*x1(i)+0.5*sqrt(3)*y1(i);
        y2(i) = 0.5*y1(i)-0.5*sqrt(3)*x1(i);
    end
end

for i = 1:numel(ts2)
    
    if atan(y3(i)/x3(i)) > pi/3
        x4(i) = 0.5*x3(i)+0.5*sqrt(3)*y3(i);
        y4(i) = 0.5*y3(i)-0.5*sqrt(3)*x3(i);
    end
end
plot3(x4,y4,ts2,'.');hold on;
plot3(x2,y2,ts,'.');
