
function Qm = Q(q,b)
    Qm = ones(q,q);
    for i = 1:q
        for j = 1:q
            if i == j
                Qm(i,j) = exp(b);
            end
        end
    end
end