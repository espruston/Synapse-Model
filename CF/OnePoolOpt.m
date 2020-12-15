x0 = [0.75,0.005,0.003,600,10]; %[p_release, k_docking, k_undocking,reserve_size, C_3]
lb = [.5, 0.0001, 0.0001,1, 10];
ub = [1, 0.01, 0.01, 1000, 10];
A = []; %p_immature - p_mature =< 0
b = [];

% opts = optimoptions('surrogateopt','CheckpointFile','checkfile2pool2.mat','PlotFcn','surrogateoptplot','UseParallel',true,'MaxTime',24*3600);
% [x,err,exitflag,output] = surrogateopt(@TwoPool2Func,lb,ub,opts);
% %[x,err,exitflag,output] = surrogateopt('checkfile2pool2.mat',opts);

problem = createOptimProblem('fmincon','objective',@OnePoolFunc,'x0',x0,'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
ms = MultiStart('FunctionTolerance',0.001,'StartPointsToRun','bounds-ineqs','UseParallel',true,'MaxTime',2*3600);
[x,err] = run(ms,problem,2);

disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), '] with an error of ', num2str(err)])