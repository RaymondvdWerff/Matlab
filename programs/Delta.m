function D = Delta(shape)
    q = shape(1);
    D = zeros(1,prod(shape));
    x = (numel(D)-1)/(q-1);
    for i = 1:q
        D(1+(i-1)*x) = 1;
    end
    D = reshape(D,shape);
end