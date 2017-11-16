format long
N = 16;

[psi0,E0] = eigs(Hamiltonian(N),1,'SA');

psi = reshape(psi0,2^(N/2),2^(N/2));
[U,s,V] = svd(psi);
pk = diag(s).^2;

psi_r = rand(2^N,1);
psi_r = psi_r/sqrt(sum(psi_r.^2));
E_r = (psi_r'*Hamiltonian(N)*psi_r)/(psi_r'*psi_r);
psi_r = reshape(psi_r,2^(N/2),2^(N/2));
[U_r,s_r,V_r] = svd(psi_r);
pk_r = diag(s_r).^2;

%{
disp('Ground state:')
disp(['The ground state energy is ',num2str(E0)])
disp(['The entanglement entropy is ',num2str(sum(-pk.*log(pk)))])

disp(' ')
disp('Random state:')
disp(['The entanglement entropy is ',num2str(sum(-pk_r.*log(pk_r)))])

semilogy(1:2^(N/2),pk)
hold on
semilogy(1:2^(N/2),pk_r)
hold off
legend('Ground state','Random state')
xlabel('number')
ylabel('pk')
axis([1 250 10^(-35) 1])
%}

Ds = [];
E0s = [];
Ers = [];
for D = 2:10:100
    Ds = [Ds,D];
    U_til = U(:,1:D);
    s_til = s(1:D,1:D);
    V_til = V(:,1:D);
    
    Ur_til = U_r(:,1:D);
    sr_til = s_r(1:D,1:D);
    Vr_til = V_r(:,1:D);
    
    psi0_til = U_til*s_til*V_til';
    psi0_til = reshape(psi0_til,2^N,1);
    E0s = [E0s,abs(E0-(psi0_til'*Hamiltonian(N)*psi0_til)/(psi0_til'*psi0_til))];
    
    psi_r_til = Ur_til*sr_til*Vr_til';
    psi_r_til = reshape(psi_r_til,2^N,1);
    Ers = [Ers,abs(E_r-(psi_r_til'*Hamiltonian(N)*psi_r_til)/(psi_r_til'*psi_r_til))];

end

semilogy(Ds,E0s)
hold on
semilogy(Ds,Ers)
hold off
legend('Ground state','Random state')
xlabel('Bond dimension')
ylabel('Energy error')


