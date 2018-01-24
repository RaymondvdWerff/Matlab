function Q = Q_Potts(q,temp,h)
    Q = ones(q)+(exp(1/temp)-1)*eye(q);
end