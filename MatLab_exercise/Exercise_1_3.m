N1 = 21;
N2 = 1000;

t1 = [0:N1-1]/(N1-1);
t2 = [0:N2-1]/(N2-1);

f1 = sin(19*pi*t1);
f2 = sin(19*pi*t2);

plot(t1,f1,'.')
hold on
plot(t2,f2)
hold off


