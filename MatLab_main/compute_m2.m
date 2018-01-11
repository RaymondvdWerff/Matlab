function [mx,my,iters,tictocs] = compute_m2(Q,q,X,tol,maxiter,ts,func)
    
    emptylist = zeros(1,numel(ts));    
    mx = emptylist;
    my = emptylist;
    iters = emptylist;
    tictocs = emptylist;
    
    delta_4D = zeros(q,q,q,q);for i=1:q; delta_4D(i,i,i,i)=1; end
    %spin1_4D = zeros(q,q,q,q);spin1_4D(1,1,1,1)=1;
    spinx_4D = zeros(q,q,q,q);for i=1:q; spinx_4D(i,i,i,i)=sin(2*pi*(i-1)/q); end
    spiny_4D = zeros(q,q,q,q);for i=1:q; spiny_4D(i,i,i,i)=cos(2*pi*(i-1)/q); end
    
    C = rand(X); C = C + C';
    T = rand(X,q,X); T = T + permute(T,[3,2,1]);
    
    for t = 1:numel(ts)
        
        temp = ts(t);
        
        Qsq = sqrtm(Q(q,temp));
        A = ncon({delta_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        %B = ncon({spin1_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        Bx = ncon({spinx_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        By = ncon({spiny_4D,Qsq,Qsq,Qsq,Qsq},{[1,2,3,4],[-1,1],[-2,2],[-3,3],[-4,4]});
        
        tic
        [C,T,iter] = func(A,C,T,X,tol,maxiter,temp);

        mx(t) = collapse(C,T,Bx)/collapse(C,T,A);
        my(t) = collapse(C,T,By)/collapse(C,T,A);
        %m(t) = sqrt(mx^2+my^2);
        %n = collapse(C,T,B)/collapse(C,T,A);
        %m(t) = (q*n-1)/(q-1);
        tictocs(t) = toc;
        iters(t) = iter;      
    end
end
    
