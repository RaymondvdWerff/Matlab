function [x2,y2] = projector(mx1,my1)
    x1 = abs(mx1);
    y1 = abs(my1);
    x2 = x1;
    y2 = y1;

    for i = 1:numel(mx1)

        if atan(y1(i)/x1(i)) > pi/3
            x2(i) = 0.5*x1(i)+0.5*sqrt(3)*y1(i);
            y2(i) = 0.5*y1(i)-0.5*sqrt(3)*x1(i);
        end
    end
end

