function Q = Q_Potts(q,temp)
    Q = ones(q)+(exp(1/temp)-1)*eye(q);
end