%Docking Increase
x0 = [0.44069,0.00079166,2.9711e-05,53.549,7.8138e-05,93.4671,13.8607]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3, C_7]
lb = [.2, 0.0002, 0.00002, 1, 0.000001, 1, 1];
ub = [.6, 0.008, 0.0002, 100, 0.0001, 100, 20];
A = [];
b = [];
fun = @Syt3DockingIncreaseFunc;

%Undocking Increase
% x0 = [0.70107,0.0030319,0.00035649,11.881,3.8596e-05,9.9347,3]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3, C_7]
% lb = [.6, 0.002, 0.0002, 2, 0.00001, 1];
% ub = [.9, 0.01, 0.002, 20, 0.0001, 10];
% A = [];
% b = [];
% fun = @Syt3UndockingFunc;

%Pool recruitment
% x0 = [0.66346,0.003346,0.00072235,7.9266,5.9651e-05,0.070014,35.809,3]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, Sty_pool_size, C_3, C_7]
% lb = [.6, 0.002, 0.00002, 1, 0.00001, 0.005, 1];
% ub = [.9, 0.01, 0.002, 20, 0.0001, 0.5, 50];
% fun = @Syt3PoolRecruitmentFunc;

% user_input_time = 12*3600;
% count = 0; count_max = 2000;
% x = x0;
% cost_best = 10;
% x_best = x0;
% tic;
% while (toc < user_input_time) && (count < count_max)
%     
%     x = lb + rand(size(x0)).*(ub-lb);
%     
%     cost = fun(x);
%     
%     if cost < cost_best
%         cost_best = cost;
%         x_best = x;
%     end
%     
%     count = count + 1;
% end


options = optimoptions('patternsearch','PlotFcn','psplotbestf','MaxFunEvals',20*3600);
[x_best,cost_best] = patternsearch(fun,x0,[],[],[],[],lb,ub,[],options);

% disp(['Best fit was [', num2str(x_best(1)), ',', num2str(x_best(2)), ',', num2str(x_best(3)), ',', num2str(x_best(4)), ',', num2str(x_best(5)), ',', num2str(x_best(6)), ',', num2str(x_best(7)), ',', num2str(x_best(8)), '] with an error of ', num2str(cost_best)])

disp(['Best fit was [', num2str(x_best(1)), ',', num2str(x_best(2)), ',', num2str(x_best(3)), ',', num2str(x_best(4)), ',', num2str(x_best(5)), ',', num2str(x_best(6)), ',', num2str(x_best(7)), '] with an cost of ', num2str(cost_best)])
