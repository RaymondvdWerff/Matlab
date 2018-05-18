function [C0,T0] = beginmatrices(Qsq,A,X,spin)
    q = size(A,1);
    
    if spin == 0
        %boundary is free
        tensor_2D = zeros(q,q);for i=1:q; tensor_2D(i,i)=1; end
        tensor_3D = zeros(q,q,q);for i=1:q; tensor_3D(i,i,i)=1; end
    else
        %boundary is fixed to given spin state
        tensor_2D = zeros(q,q);tensor_2D(spin,spin) = 1;
        tensor_3D = zeros(q,q,q);tensor_3D(spin,spin,spin)=1;
    end
    
    C0 = ncon({tensor_2D,Qsq,Qsq},{[1,2],[-1,1],[-2,2]});
    T0 = ncon({tensor_3D,Qsq,Qsq,Qsq},{[1,2,3],[-1,1],[-2,2],[-3,3]});
    
    while size(T0,1) < X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TA = ncon({T0,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TA},{[-1,1,2],[1,2,-2,-3,-4]});
        C0 = reshape(M,q*size(T0,1),q*size(T0,1));
        T0 = reshape(TA,[q*size(T0,1),q,q*size(T0,1)]);
        C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
        T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
    end
    
%     if size(T0,1) > X
%         [U,~,~] = svd(C0);
%         U_til = U(:,1:X);
%         C0 = ncon({C0,U_til,U_til},{[1,2],[1,-1],[2,-2]});
%         T0 = ncon({T0,U_til,U_til},{[1,-2,2],[1,-1],[2,-3]});
%         C0 = (C0+permute(C0,[2,1]))./max(abs(C0(:)));
%         T0 = (T0+permute(T0,[3,2,1]))./max(abs(T0(:)));
%     end
    
    for i = 1:X
        CT = ncon({C0,T0},{[1,-2],[-1,-3,1]});
        TA = ncon({T0,A},{[-1,1,-4],[1,-2,-3,-5]});
        M = ncon({CT,TA},{[-1,1,2],[1,2,-2,-3,-4]});
        
        [U,~,~] = tensorsvd(M,[1,2],[3,4],X);
        
        C0 = ncon({U,M,U},{[1,2,-1],[1,2,3,4],[3,4,-2]});
        C0 = (C0 + permute(C0,[2,1]))/max(C0(:));
        T0 = ncon({U,TA,U},{[1,2,-1],[1,2,-2,3,4],[3,4,-3]});
        T0 = (T0 + permute(T0,[3,2,1]))/max(T0(:));
    end
end