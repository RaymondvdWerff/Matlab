x = 0;
y = 0;
dt = 0.1;

for i = 0:1000
    vx = randn(1);
    vy = randn(1);
    x = x + vx*dt;
    y = y + vy*dt;
    plot(x,y,'.')
    hold on
    %pause(0.01)
end