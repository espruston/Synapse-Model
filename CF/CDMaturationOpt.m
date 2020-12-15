x0 =[0.95,0,0.0019237,0.0006432,0.0037323,0.0077772,10];
%initial guess p_mature = x(1); p_immature = x(2); k_docking = x(3); k_undocking = x(4); k_maturation = x(5); k_dematuration = x(6);
A = [-1 1 0 0 0 0 0]; %p_immature - p_mature =< 0
b = 0;
Aeq = [];
beq = [];
lb = [.6, .0, .0001, 0.0001, .0001, 0.0001, 1];
ub = [1, .2, .01, .01, .01, .01, 10];
nonlcon = [];


opts = optimoptions('surrogateopt','CheckpointFile','checkfileCDMat.mat','PlotFcn','surrogateoptplot','UseParallel',true,'MaxTime',3600);
[x,err,exitflag,output] = surrogateopt(@CDMaturationFunc,lb,ub,opts);
%[x,err,exitflag,output] = surrogateopt('checkfileCDMat.mat',opts);


% problem = createOptimProblem('fmincon','objective',@CDMaturationFunc,'x0',x0,'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
% ms = MultiStart('FunctionTolerance',0.001,'StartPointsToRun','bounds-ineqs','UseParallel',true,'MaxTime',7*3600);
% [x,err] = run(ms,problem,7);


disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), ',', num2str(x(7)), '] with an error of ', num2str(err)])