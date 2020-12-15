x0 = [0.8,0.001,0.0001,50,0.0002,5e6]; %[p_release, k_docking, k_undocking,reserve_size, k_refill, C_Ca]
lb = [.5, 0.0001, 0.0001, 1, 0.0001, 1e6];
ub = [1, 0.01, 0.01, 100, 0.001, 1e7];
A = [];
b = [];

opts = optimoptions('surrogateopt','CheckpointFile','checkfileRefill.mat','PlotFcn','surrogateoptplot','UseParallel',true,'MaxTime',24*3600);
[x,err,exitflag,output] = surrogateopt(@ReserveRefillFunc,lb,ub,opts);
%[x,err,exitflag,output] = surrogateopt('checkfileRefill.mat',opts);

% problem = createOptimProblem('fmincon','objective',@ReserveRefillFunc,'x0',x0,'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
% ms = MultiStart('FunctionTolerance',0.1,'StartPointsToRun','bounds-ineqs','UseParallel',true,'MaxTime',600);
% [x,err] = run(ms,problem,1);

disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), '] with an error of ', num2str(err)])
