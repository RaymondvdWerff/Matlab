function A = get_symtensor(c)
    % create a D=2 rotational symmetric tensor
    % depending on 12 parameters c
    A=zeros(2,2,2,2,2);
    a=1;b=2;
    A(a,a,a,a,a) = c(1);
    A(a,b,a,a,a) = c(2);
    A(a,a,b,a,a) = c(2);
    A(a,a,a,b,a) = c(2);
    A(a,a,a,a,b) = c(2);
    A(a,b,b,a,a) = c(3);
    A(a,a,b,b,a) = c(3);
    A(a,a,a,b,b) = c(3);
    A(a,b,a,a,b) = c(3);
    A(a,b,a,b,a) = c(4);
    A(a,a,b,a,b) = c(4);
    A(a,b,b,b,a) = c(5);
    A(a,a,b,b,b) = c(5);
    A(a,b,a,b,b) = c(5);
    A(a,b,b,a,b) = c(5);
    A(a,b,b,b,b) = c(6);
    
    a=2;b=1;
    A(a,a,a,a,a) = c(7);
    A(a,b,a,a,a) = c(8);
    A(a,a,b,a,a) = c(8);
    A(a,a,a,b,a) = c(8);
    A(a,a,a,a,b) = c(8);
    A(a,b,b,a,a) = c(9);
    A(a,a,b,b,a) = c(9);
    A(a,a,a,b,b) = c(9);
    A(a,b,a,a,b) = c(9);
    A(a,b,a,b,a) = c(10);
    A(a,a,b,a,b) = c(10);
    A(a,b,b,b,a) = c(11);
    A(a,a,b,b,b) = c(11);
    A(a,b,a,b,b) = c(11);
    A(a,b,b,a,b) = c(11);
    A(a,b,b,b,b) = c(12);
end