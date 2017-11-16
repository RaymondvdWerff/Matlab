A = rand(3,2,5,4);
B = rand(2,5,3);

A = permute(A,[1,2,4,3]);
B = permute(B,[2,1,3]);

sA = size(A);
Adimleft = sA(1:3);
Adimright = sA(4);

sB = size(B);
Bdimleft = sB(1);
Bdimright = sB(2:3);

MA = reshape(A,prod(Adimleft),Adimright);
MB = reshape(B,Bdimleft,prod(Bdimright));

MC = MA*MB;
C = reshape(MC,[Adimleft,Bdimright]);

size(C)