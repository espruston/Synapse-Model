x0 = [0.72436,0.0023773,0.0005734,8.9984,0.0001,100000.0001]; %[p_release, k_docking, k_undocking,reserve_size, k_refill, C_Ca]
lb = [.5, 0.001, 0.0001, 1, 0.0001, 1e5];
ub = [1, 0.01, 0.001, 20, 0.0001, 1e8];
A = [];
b = [];

% opts = optimoptions('surrogateopt','CheckpointFile','checkfileRefill.mat','PlotFcn','surrogateoptplot','UseParallel',false,'MaxTime',3600);
% [x,err,exitflag,output] = surrogateopt(@ReserveRefillFunc,lb,ub,opts);
%[x,err,exitflag,output] = surrogateopt('checkfileRefill.mat',opts);

% problem = createOptimProblem('fmincon','objective',@ReserveRefillFunc,'x0',x0,'lb',lb,'ub',ub,'Aineq',A,'bineq',b);
% ms = MultiStart('FunctionTolerance',0.005,'StartPointsToRun','bounds-ineqs','UseParallel',false,'MaxTime',24*3600);
% [x,err] = run(ms,problem,10);

%random search
% range = ub-lb;
% user_input_time = 54000;
% count = 0; count_max = 5000;
% x_best = [];
% err_best = inf;
% tic;
% while (toc < user_input_time) && (count < count_max) % always include a failsafe!
%     count = count + 1;
%     x = lb+rand(1,6).*range;
%     err = ReserveRefillFunc(x);
%     
%     if err < err_best
%         err_best = err;
%         x_best = x;
%     end
%     
% end
% disp([num2str(count),' parameter sets tested']);

% fun = @ReserveRefillFunc;
% options = optimset('PlotFcn','optimplotfval','MaxFunEvals',1000,'Display','iter');
% [x_best,err_best,eflag,output] = fminunc(fun,x0,options);

fun = @ReserveRefillFunc;
options = optimoptions('patternsearch','PlotFcn','psplotbestf','MaxFunEvals',500);
[x_best,err_best] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);

disp(['Best fit was [', num2str(x_best(1)), ',', num2str(x_best(2)), ',', num2str(x_best(3)), ',', num2str(x_best(4)), ',', num2str(x_best(5)), ',', num2str(x_best(6)), '] with an error of ', num2str(err_best)])

%disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), '] with an error of ', num2str(err)])
