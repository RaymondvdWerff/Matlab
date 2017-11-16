x = 0;
y = 0;
vx = 0.1;
vy = 0.2;
dt = 0.1;

for i = 1:500
    if (x + vx*dt < -1)|(x + vx*dt >1)
        vx = -vx;
    end
    if (y + vy*dt < -1)|(y + vy*dt > 1)
        vy = -vy;
    end
    x = x + vx*dt;
    y = y + vy*dt;
    
    plot([-1 1],[1 1],'b-')
    hold on
    axis([-1.1 1.1 -1.1 1.1])
    plot([-1 1],[-1 -1],'b-')
    plot([1 1],[-1 1],'b-')
    plot([-1 -1],[-1 1],'b-')
    plot(x,y,'*');
    hold off
    pause(0.01);
end


