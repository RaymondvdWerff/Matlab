function [Al,lambda,counter] = MPOFixedPoint2(T,Ar,m,eta)

    [Al,C] = LeftOrthonormalize2(Ar,eye(size(Ar,2)),m,eta);
    delta = 10;
    maxiter = 1000;
    counter = 0;
    
    while delta > eta
    
        Ml = ncon({Al,conj(Al),T},{[1,-4,-1],[2,-6,-3],[1,2,-5,-2]});
        siMl = size(Ml);
        Ml = reshape(Ml,prod(siMl(1:3)),prod(siMl(4:6)));
        Mr = ncon({Ar,conj(Ar),T},{[1,-1,-4],[2,-3,-6],[1,2,-2,-5]});
        siMr = size(Mr);
        Mr = reshape(Mr,prod(siMr(1:3)),prod(siMr(4:6)));

        [lambda,Fl] = LeadingEigenvector(Ml,rand(size(Ml,2),1),m,delta/10);
        [~,Fr] = LeadingEigenvector(Mr,rand(size(Mr,2),1),m,delta/10);
        Fl = reshape(Fl,[siMl(4),siMl(5),siMl(6)]);
        Fr = reshape(Fr,[siMr(4),siMr(5),siMr(6)]);
        
        Fl =  Fl./ncon({Fl,Fr,C,conj(C)},{[1,3,2],[4,3,5],[1,4],[2,5]});
    
        Ac = ncon({C,Ar},{[-2,1],[-1,1,-3]});
        
        M = ncon({Fl,Fr,T},{[-5,1,-2],[-6,2,-3],[-4,-1,1,2]});
        siM = size(M);
        M = reshape(M,prod(siM(1:3)),prod(siM(4:6)));
        [~,Ac] = LeadingEigenvector(M,reshape(Ac,size(M,2),1),m,delta/10);
        Ac = reshape(Ac,[siM(4),siM(5),siM(6)]);
        
        M = ncon({Fl,Fr},{[-3,1,-1],[-4,1,-2]});
        siM = size(M);
        M = reshape(M,prod(siM(1:2)),prod(siM(3:4)));
        [~,C] = LeadingEigenvector(M,reshape(C,size(M,2),1),m,delta/10);
        C = reshape(C,siM(3),siM(4));
        
        siAc = size(Ac);
        [Qa,Ra] = qrpos(reshape(Ac,prod(siAc(1:2)),siAc(3)));
        Qa = reshape(Qa,[siAc(1),siAc(2),siAc(3)]); 
        [Qc,Rc] = qrpos(C);
        Al = ncon({Qa,conj(Qc)},{[-1,-2,1],[-3,1]});
        
        Aclq = permute(Ac,[2,3,1]);
        siAclq = size(Aclq);
        [La,Qa] = lqpos(reshape(Aclq,siAclq(1),prod(siAclq(2:3))));
        Qa = reshape(Qa,[siAclq(1),siAclq(2),siAclq(3)]);
        Qa = permute(Qa,[3,1,2]);
        [Lc,Qc] = lqpos(C);
        Ar = ncon({conj(Qc),Qa},{[1,-2],[-1,1,-3]});
        
        delta = max(sqrt(ncon({Ra-Rc,conj(Ra-Rc)},{[1,2],[1,2]})),sqrt(ncon({La-Lc,conj(La-Lc)},{[1,2],[1,2]})))
        
        counter = counter + 1;

        if counter > maxiter
            disp('Not converged');
            break;
        end

    end
end