N = 1000;
t_begin = 0;
t_end = 1;
t_step = (t_end-t_begin)/N;

t = [t_begin:t_step:t_end];
f1 = Sinat(2,t);
f2 = Sinat(4,t);
f3 = Sinat(8,t);

g1 = zeros(N+1);
g2 = zeros(N+1);
g3 = zeros(N+1);

for i = 1:N
    g1(i) = sum((f1(1:i)+f1(2:i+1)))*t_step/2;
    g2(i) = sum((f2(1:i)+f2(2:i+1)))*t_step/2;
    g3(i) = sum((f3(1:i)+f3(2:i+1)))*t_step/2;
end

plot(t,g1)
hold on
plot(t,g2)
plot(t,g3)
hold off
legend('sin(2*pi*t)','sin(4*pi*t)','sin(8*pi*t)')
    