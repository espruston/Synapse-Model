x0 = [0.81832,0.01246,0.0019237,0.0016432,0.0042675,0.0037772,9.2913]; %initial guess p_mature = x(1); p_immature = x(2); k_docking = x(3); k_undocking = x(4); k_maturation = x(5); k_dematuration = x(6);
A = [-1 1 0 0 0 0 0]; %p_immature - p_mature =< 0
b = 0;
Aeq = [];
beq = [];
lb = [.6, .0, .0001, 0.0001, .0001, 0.0001, 1];
ub = [1, .2, .01, .01, .01, .01, 10];
nonlcon = [];


%Local min solver
% options = optimoptions('patternsearch','UseParallel',true,'MaxTime',7200);
% [x,err] = patternsearch(@MaturationCDRFunc,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);


%options = optimoptions('ga','UseParallel',true,'MaxTime',104400);
%[x,err] = ga(@MaturationCDRFunc,7,A,b,Aeq,beq,lb,ub,nonlcon,options);


% opts = optimoptions('surrogateopt','CheckpointFile','checkfile.mat','PlotFcn','surrogateoptplot','UseParallel',true,'MaxTime',3600+54000);
% [x,err,exitflag,output] = surrogateopt(@MaturationCDRFunc,lb,ub,opts);
% [x,err,exitflag,output] = surrogateopt('checkfile.mat',opts);



problem = createOptimProblem('fmincon','objective',@MaturationCDRFunc,'x0',x0,'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
ms = MultiStart('FunctionTolerance',0.001,'StartPointsToRun','bounds-ineqs','UseParallel',true,'MaxTime',7*3600);
[x,err] = run(ms,problem,7);


disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), ',', num2str(x(7)), '] with an error of ', num2str(err)])