t = [0:1000]/1000;
y1 = sin(2*pi*t);
y2 = exp(-5*t);
y3 = cos(t).*exp(-t);

plot(t,y1)
hold on
plot(t,y2)
plot(t,y3)
hold off
xlabel('t')
ylabel('y')
legend('sin(2*pi*t)','exp(-5*t)','cos(t)*exp(-t)') 