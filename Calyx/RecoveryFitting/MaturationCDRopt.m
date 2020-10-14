x0 = [0.51061,0.096496,0.010349,0.00079271,0.00046181,0.00010257,5.3567]; %initial guess p_mature = x(1); p_immature = x(2); k_docking = x(3); k_undocking = x(4); k_maturation = x(5); k_dematuration = x(6); C_3 = x(7)
A = [-1 1 0 0 0 0 0; 0 0 -1 1 0 0 0; 0 0 0 0 -1 1 0]; %p_immature - p_mature =< 0, k_undocking - k_docking =< 0, k_dematuration - k_maturation =< 0  
b = [0; 0; 0];
Aeq = [];
beq = [];
lb = [.1, .01, .0001, 0, .00001, 0, 0];
ub = [.7, .2, .01, .001, .01, .001, 10];
nonlcon = [];
%options = optimoptions('patternsearch','UseParallel',true,'MaxTime',60000);
opts = optimoptions('surrogateopt','CheckpointFile','checkfile.mat','UseParallel',true,'MaxTime',10800+10800,'MaxFunctionEvaluations',5000);
%opts.MaxTime = 1200;
%opts.UseParallel = true;

%[x,err] = patternsearch(@MaturationCDRFunc,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
[x,err,exitflag,output] = surrogateopt('checkfile.mat',opts);
%[x,err,exitflag,output] = surrogateopt(@MaturationCDRFunc,lb,ub,opts);

disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), ',', num2str(x(7)), '] with an error of ', num2str(err)])