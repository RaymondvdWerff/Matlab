function y = m_exact(x)
    if x < 1/log(1+sqrt(2))
        y = (1-sinh(1/x)^(-4))^(1/8);
    else
        y = 0;
    end
end
