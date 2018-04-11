function Q = Q_trIsing(q,temp,h)
    b = 1/temp;
    Q = ones(4);
    Q(1,1) = exp(b);Q(2,2) = Q(1,1);
    Q(1,2) = exp(-b);Q(2,1) = Q(1,2);
    Q(1,3) = exp(b*h/4);Q(3,1) = Q(1,3);Q(2,3) = Q(1,3);Q(3,2) = Q(1,3);
    Q(1,4) = exp(-b*h/4);Q(4,1) = Q(1,4);Q(2,4) = Q(1,4);Q(4,2) = Q(1,4);
    Q(3,3) = exp(b*h/2);
    Q(4,4) = exp(-b*h/2);
end