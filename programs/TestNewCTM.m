q = 6;
X = 10;
tol = 1e-4;
maxiter = 10;

temp = 0.8;

delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
Qsq = sqrtm(Q_clock(q,temp,0));
A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
    


for iter = 1:maxiter
    
    M = ncon({A,A,A,A},{[-1,2,-5,1],[-2,4,1,-7],[2,-3,-6,3],[4,-4,3,-8]});
    
    [U,s,~] = tensorsvd(M,[1,2],[3,4,5,6,7,8],X);
    s = s/max(s(:));
    
    A = ncon({M,U,U,U,U},{[1,2,3,4,5,6,7,8],[1,2,-1],[3,4,-2],[5,6,-3],[7,8,-4]});
    A = A/max(A(:));

    if iter > 1
        delta = sum(sum(abs(s-sold)))
        if delta < tol
            break;
        end
    end  
    sold = s;           
end