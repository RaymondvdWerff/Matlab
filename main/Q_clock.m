function Q = Q_clock(q,temp,h)
    Q = zeros(q);
    for i = 1:q
        for j = 1:q
            Q(i,j) = exp((cos(2*pi*abs(i-j)/q)+h*(cos(2*pi*(i-1)/q)+cos(2*pi*(j-1)/q))/4)/temp);
        end
    end
end