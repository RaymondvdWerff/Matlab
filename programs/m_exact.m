function y = m_exact(ts)
    
    y = zeros(1,numel(ts));
    
    for temp = 1:numel(ts)
        t = ts(temp);
        if t < 2/log(1+sqrt(2))
            y(temp) = (1-sinh(2/t)^(-4))^(1/8);
        else
            y(temp) = 0;
        end
    end
end
