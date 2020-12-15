x0 = [0.6,0.005,0.003,0.001,0.003,0.75,10]; %[p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, size_rel, C_3]
lb = [.5, 0.0001, 0.0001, 0.0001, 0.0001, 0.25, 10];
ub = [1, 0.01, 0.01, 0.01, 0.01, 1, 10];
A = [-1 1 0 0 0 0 0]; %p_immature - p_mature =< 0
b = 0;

% opts = optimoptions('surrogateopt','CheckpointFile','checkfile2pool2.mat','PlotFcn','surrogateoptplot','UseParallel',true,'MaxTime',24*3600);
% [x,err,exitflag,output] = surrogateopt(@TwoPool2Func,lb,ub,opts);
% %[x,err,exitflag,output] = surrogateopt('checkfile2pool2.mat',opts);

problem = createOptimProblem('fmincon','objective',@TwoPool2Func,'x0',x0,'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
ms = MultiStart('FunctionTolerance',0.001,'StartPointsToRun','bounds-ineqs','UseParallel',true,'MaxTime',12*3600);
[x,err] = run(ms,problem,12);

disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), ',', num2str(x(7)),'] with an error of ', num2str(err)])