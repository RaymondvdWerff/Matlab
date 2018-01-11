N = 1000;
t = [0:N]/N;
f = zeros(N+1,1);

for i = 1:N
    if exp(5*t(i))<10
        f(i) = exp(5*t(i));
    else
        f(i) = 10;
    end
end
    
plot(t,f)
axis([0 1 0 11])
xlabel('t')
ylabel('f')