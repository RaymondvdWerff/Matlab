
%[Al,C1,counter1] = LeftOrthonormalize2(A,eye(size(A,2)),eta);

%[Al,lambda,counter] = MPOFixedPoint3(T,A,m,eta);
%counter


%[Al,C1,counter1] = LeftOrthonormalize1(A,eye(size(A,2)),eta);
%[Ar,C2,counter2] = LeftOrthonormalize1(permute(A,[1,3,2]),eye(size(A,2)),eta);

A = rand(2,3,3);
A = A + permute(A,[1,3,2]);
eta = 10^(-10);

[Al,C1] = LeftOrthonormalize1(A,eye(size(A,2)),eta);
[Ar,C2] = LeftOrthonormalize1(permute(A,[1,3,2]),eye(size(A,2)),eta);

Ar = permute(Ar,[1,3,2]);
C2 = permute(C2,[2,1]);

[Al_2,C1_2] = LeftOrthonormalize3(A,eye(size(A,2)),eta);
[Ar_2,C2_2] = LeftOrthonormalize3(permute(A,[1,3,2]),eye(size(A,2)),eta);

Ar_2 = permute(Ar_2,[1,3,2]);
C2_2 = permute(C2_2,[2,1]);

Ac_2 = ncon({C1_2,A,C2_2},{[-2,1],[-1,1,2],[2,-3]});
Ac_2 = Ac_2./max(abs(Ac_2(:)));

%{
Ac1 = ncon({Al,C1,C2},{[-1,-2,1],[1,2],[2,-3]});
Ac2 = ncon({C1,C2,Ar},{[-2,1],[1,2],[-1,2,-3]});
Ac1 = Ac1/max(abs(Ac1(:)));
Ac2 = Ac2/max(abs(Ac2(:)));

delta1 = sqrt(ncon({Ac1-Ac2,conj(Ac1-Ac2)},{[1,2,3],[1,2,3]})) 
delta2 = sqrt(ncon({Ac-Ac1,conj(Ac-Ac1)},{[1,2,3],[1,2,3]})) 
delta3 = sqrt(ncon({Ac-Ac2,conj(Ac-Ac2)},{[1,2,3],[1,2,3]})) 
%}

%{
A1 = ncon({C1,A},{[-2,1],[-1,1,-3]});
A2 = ncon({Al,C1},{[-1,-2,1],[1,-3]});
A1 = A1/max(abs(A1(:)));
A2 = A2/max(abs(A2(:)));
delta = sqrt(ncon({A1-A2,conj(A1-A2)},{[1,2,3],[1,2,3]})) 
%} 

%{
R = ncon({conj(C2),C2},{[-2,1],[-1,1]});
R = R./max(abs(R(:)))

R_tilde = ncon({R,A,conj(A)},{[1,2],[3,-1,1],[3,-2,2]});
R_tilde = R_tilde./max(abs(R_tilde(:)))

L = ncon({conj(C1),C1},{[1,-2],[1,-1]});
L = L./max(abs(L(:)))

L_tilde = ncon({L,A,conj(A)},{[1,2],[3,1,-1],[3,2,-2]});
L_tilde = L_tilde./max(abs(L_tilde(:)))
%}

%{
Oz = [1 0;0 -1];
exp1 = ncon({C1,C2,Ar,conj(C1),conj(C2),conj(Ar),Oz},{[2,3],[3,4],[5,4,7],[2,1],[1,8],[6,8,7],[5,6]});
nrm1 = ncon({C1,C2,Ar,conj(C1),conj(C2),conj(Ar)},{[2,3],[3,4],[5,4,6],[2,1],[1,7],[5,7,6]});
exp1/nrm1

exp2 = ncon({C1,C2,Al,conj(C1),conj(C2),conj(Al),Oz},{[8,1],[1,2],[6,7,8],[4,3],[3,2],[5,7,4],[6,5]});
nrm2 = ncon({C1,C2,Al,conj(C1),conj(C2),conj(Al)},{[7,1],[1,2],[5,6,7],[4,3],[3,2],[5,6,4]});
exp2/nrm2

exp3 = ncon({Ac,Oz,conj(Ac)},{[1,2,4],[1,3],[3,2,4]});
nrm3 = ncon({Ac,conj(Ac)},{[1,2,3],[1,2,3]});
exp3/nrm3

exp4 = ncon({C1,A,A,C2,conj(C1),conj(A),conj(A),conj(C2),Oz},{[1,3],[5,3,7],[9,7,10],[10,2],[1,4],[6,4,8],[9,8,11],[11,2],[5,6]});
nrm4 = ncon({C1,A,A,C2,conj(C1),conj(A),conj(A),conj(C2)},{[1,3],[5,3,6],[8,6,9],[9,2],[1,4],[5,4,7],[8,7,10],[10,2]});
exp4/nrm4
%}



