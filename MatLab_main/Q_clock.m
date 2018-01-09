function Q = Q_clock(q,temp)
    Q = zeros(q);
    for i = 1:q
        for j = 1:q
            Q(i,j) = exp(cos(2*pi*abs(i-j)/q)/temp);
        end
    end
end