function [lambda,vector] = LeadingEigenvector(A,b,m,tol)

    H = zeros(m+1,m);
    V = zeros(size(A,2),m);
    V(:,1) = b/norm(b);
    
    for j = 1:m
        w = A*V(:,j);
        for i = 1:j
            H(i,j) = w'*V(:,i);
            w = w - H(i,j)*V(:,i);
        end
        
        H(j+1,j) = norm(w);
        
        if H(j+1,j) < tol
            break
        end
        
        if j < m
            V(:,j+1) = w./H(j+1,j);
        end
    end
    
    H = H(1:j,1:j);
    V = V(:,1:j);
    [vector,lambda] = eigs(H,1,'LM');
    vector = V*vector;

end