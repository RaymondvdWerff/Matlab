q = 6;
tol = 6;
X = 20:5:200;
maxiter = 2000;

t_begin = 0.88;
t_step = 0.001;
t_end = 0.94;

ts = t_begin:t_step:t_end;
xs = zeros(numel(X),numel(ts));
ys = zeros(numel(X),numel(ts));
iters = zeros(numel(X),numel(ts));
tictocs = zeros(numel(X),numel(ts));

for x = 1:numel(X)
    disp(['X = ' , num2str(X(x))]);
    [ev2,ev3,iter,tictoc] = compute_corr2(@Q_clock,q,X(x),10^(-tol),maxiter,ts,@FPCM);
    xs(x,:) = log(ev2/ev3);
    ys(x,:) = -log(ev2);
    iters(x,:) = iter;
    tictocs(x,:) = tictoc;
end

save('C_corr_q6X20-200tol6_FPCM.mat_r','X','ts','xs','ys','iters','tictocs');