function [HL,HM,HR] = HMPO(h)
    
    Sz = [1 0;0 -1];
    Sx = [0 1;1 0];
    
    HL = zeros(2,2,3);
    HL(:,:,1) = -h*Sx;
    HL(:,:,2) = Sz;
    HL(:,:,3) = eye(2);
    
    HM = zeros(2,2,3,3);
    HM(:,:,1,1) = eye(2);
    HM(:,:,2,1) = -Sz;
    HM(:,:,3,1) = -h*Sx;
    HM(:,:,3,2) = Sz;
    HM(:,:,3,3) = eye(2);
    
    HR = zeros(2,2,3);
    HR(:,:,1) = eye(2);
    HR(:,:,2) = -Sz;
    HR(:,:,3) = -h*Sx;
    
end
    