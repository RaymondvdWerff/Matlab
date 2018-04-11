q = 2;
tol = 6;
X = 20:5:100;
maxiter = 20000;

t_begin = 2.267;
t_step = 0.0001;
t_end = 2.273;

ts = t_begin:t_step:t_end;
xs = zeros(numel(X),numel(ts));
ys = zeros(numel(X),numel(ts));
iters = zeros(numel(X),numel(ts));
tictocs = zeros(numel(X),numel(ts));

for x = 1:numel(X)
    disp(['X = ' , num2str(X(x))]);
    [ev2,ev3,iters1,tictocs1] = compute_corr2(@Q_clock,q,X(x),10^(-tol),maxiter,ts,@CTM);
    xs(x,:) = log(ev2/ev3);
    ys(x,:) = -log(ev2);
    iters(x,:) = iters1;
    tictocs(x,:) = tictocs1;
end

save('C_corr_q2X20-100tol6_FPCM.mat','X','ts','xs','ys','iters','tictocs');