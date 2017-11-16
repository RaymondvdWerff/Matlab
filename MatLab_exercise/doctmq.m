function E=doctmq(A, H, chi, tol, verbose)
%Use corner-transfer matrix method to compute environment of a 2D quantum state given by tensor A
%(rotational and translational symmetry is assumed). H is the Hamiltonian

if nargin<1 
    %default values for testing:
    lambda=3;
    H = get_H_trans_ising(lambda);
    % get a symmetric tensor parameterized by 12 values
    A = get_symtensor([1 0.5 0.3 0.5 0.1 0.1 1 0.2 0.3 0.5 0.1 1]);
    chi = 8;
    tol = 1e-7;
end
if nargin<4
    verbose=1;
end

%more parameters:
maxsteps=20000;
minsteps=4;
seed=1;

D=size(A,2);
D2=D*D;

%initialize symmetric random boundary
rand('seed',seed);
chi0=2;
C=rand(chi0,chi0);
C=C+permute(C,[2 1]);
T=rand(chi0,chi0,D,D);     %first two legs are boundary indices
T=T+permute(T,[2 1 3 4]);
T=T+permute(T,[1 2 4 3]);
T=reshape(T,chi0,chi0,D2); %fuse legs together

% another possibility: initialize boundary from A tensor
% C=ncon({A,conj(A)},{[1 2 3 -1 -3],[1 2 3 -2 -4]});
% C=reshape(C,[D2,D2]);
% T=ncon({A,conj(A)},{[1 2 -5 -1 -3],[1 2 -6 -2 -4]});
% T=reshape(T,[D2,D2,D2]);

%compute reduced tensor
aa = ncon({A,conj(A)},{[1 -1 -3 -5 -7],[1 -2 -4 -6 -8]});
aa = reshape(aa, [D2,D2,D2,D2]);

%fprintf('\nDo contraction for beta=%g, chi=%g, tol=%g, minsteps=%g, maxsteps=%g \n',beta,chi,tol,minsteps,maxsteps);
sold = zeros(chi0,1);

% do main loop
for i=1:maxsteps
    % do a ctm iteration and get new corner spectrum
    snew = doctmstep();
    %compute difference in corner spectrum
    if numel(sold) == numel(snew)
        diffs =  norm(snew-sold);
    else
        diffs = inf;
    end
    
    if verbose
        %compute energy
        E = getenergy();
        fprintf('step:%d,  E=%g,  diffs=%g \n',i,E,diffs);
    end
    if (diffs<tol)  && (i>minsteps) 
        break;
    end
    sold = snew;
end


%final calculation:
E = getenergy();
fprintf('chi=%d,  E=%g\n',chi, E);

function s=doctmstep()
    % do one CTM step and return the new corner spectrum
    
    % compute corner (e.g. upper left)
    nC=ncon({C,T,T,aa},{[1 2],[-1 1 3],[2 -3 4],[-2 3 4 -4]});
    
    % compute isometry
    si=size(nC);
    % reshape tensor into a matrix
    mC=reshape(nC,prod(si(1:2)),prod(si(3:4)));
    %symmetrize (T might be non-symmetric due to round off errors)
    mC=(mC+mC')/2;
    cchi=min(chi,size(mC,1));
    
    %get chi largest eigenstates
%     [U,s]=eig(mC);
%     [sd,ii]=sort(abs(diag(s)),'descend');
%     s=s(ii,ii);U=U(:,ii);
    [U,s,Vt]=svd(mC,'econ');
    U=U(:,1:cchi);
    U=reshape(U,[si(1) si(2) cchi]);
    s=diag(s);
    s=s(1:cchi);
    s=s/max(s);
    
    % Compute renormalized corner
    % note that since everything is real conj(U)=U
    C=ncon({nC,conj(U),U},{[1 2 3 4],[1 2 -1 ],[3 4 -2]});
    C=normalizetensor(C);
    % the result is equal to s here: C = diag(s);
    
    % new edge
    T=ncon({U,T,aa,conj(U)},{[1 2 -1],[1 4 3],[-3 2 3 5],[4 5 -2]});  
    T=normalizetensor(T);
end

function E=getenergy()
    %compute 2-site reduced density matrix and compute energy
    %split D2 index of edge tensor (to connect to ket and bra tensor)
    s = size(T);
    T2 = reshape(T,[s(1),s(2),D,D]);
    
    leftpart = ncon({C,T2,C,T2,A,conj(A),T2},{[3 1],[1 2 4 6],[10 2],[-1 3 5 7],[-5 -2 5 4 8],[-6 -3 7 6 9],[-4 10 8 9]}); 
    rightpart = leftpart;
    rho2 = ncon({leftpart,rightpart},{[1 2 3 4 -1 -3],[1 2 3 4 -2 -4]});
    
    nrm = ncon({rho2},{[1 2 1 2]});
    E = ncon({rho2,H},{[3 4 1 2],[1 2 3 4]});
    E = E/nrm;
end
end

function T=normalizetensor(T)
    %normalize a tensor by dividing by its largest number
    T=T/max(abs(T(1:end)));
end