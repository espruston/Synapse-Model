%Docking Increase
% x0 = [0.69163,0.0029594,0.00021448,13.2628,3.0985e-05,9.5945]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3]
% lb = [.6, 0.002, 0.0002, 2, 0.00001, 1];
% ub = [.9, 0.01, 0.002, 20, 0.0001, 10];
% A = [];
% b = [];
% fun = @Syt3DockingIncreaseFunc;

%Docking and undocking Increase
% x_best = [0.737502,0.002,0.00146176,14.2557,3.5552e-05,6.48438]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3]
% lb = [.4, 0.0002, 0.0002, 1, 0.00001, 1];
% ub = [1, 0.9, 0.9, 20, 0.9, 50];
% A = [];
% b = [];
% fun = @Syt3DockingAndUndockingIncreaseFunc;

%Docking Increase 2 State
% x0 = [0.671166,0.599586,0.00569202,0.00572867,0.137335,0.0176688,27.407,0.0832238,30.8695]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3]
% lb = [0, 0, 0.0001, 0.0001, .0001, .00001, 10, 0.00001, 0];
% ub = [1, 1, 0.1, 0.1, 0.9, 0.9, 50, 0.1, 50];
% A = [];
% b = [];
% fun = @DockingIncrease2StateFunc;

%Undocking Increase
% x0 = [0.70107,0.0030319,0.00035649,11.881,3.8596e-05,9.9347]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3]
% lb = [.6, 0.002, 0.0002, 2, 0.00001, 1];
% ub = [.9, 0.01, 0.002, 20, 0.0001, 10];
% A = [];
% b = [];
% fun = @Syt3UndockingFunc;

%Pool recruitment
% x0 = [0.66346,0.003346,0.00072235,7.9266,5.9651e-05,0.070014,35.809]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, Sty_pool_size, C_3]
% lb = [.6, 0.002, 0.00002, 1, 0.00001, 0.005, 1];
% ub = [.9, 0.01, 0.002, 20, 0.0001, 0.5, 50];
% fun = @Syt3PoolRecruitmentFunc;

%P increase
x_best = [0.75789,0.0037549,0.0021137,3.0745,0.00060017,1.89]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3]
lb = [0, 0.0002, 0.0002, 1, 0.00001, 1];
ub = [1, 0.9, 0.9, 20, 0.9, 50];
A = [];
b = [];
fun = @Syt3PIncreaseFunc;

%Basic
% x_best = [0.688,0.0030421,0.0032928,9.5651,6.607e-05]; %[p_release, k_docking, k_undocking, reserve_size, k_refill]
% lb = [.4, 0.0002, 0.0002, 1, 0.00001];
% ub = [1, 0.09, 0.09, 20, 0.9];
% A = [];
% b = [];
% fun = @WTBasicFunc;

% user_input_time = 1*3600;
% count = 0; count_max = 2000;
% x = x_best;
% cost_best = 10;
% tic;
% while (toc < user_input_time) && (count < count_max)
%     
%     x = lb + rand(size(lb)).*(ub-lb);
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


options = optimoptions('patternsearch','PlotFcn','psplotbestf','MaxIterations',1e6,'MaxTime',.15*3600);
[x_best,cost_best] = patternsearch(fun,x_best,[],[],[],[],lb,ub,[],options);

disp(['Best fit was x = [', num2str(x_best), '] with an cost of ', num2str(cost_best)])

% disp(['Best fit was [', num2str(x_best(1)), ',', num2str(x_best(2)), ',', num2str(x_best(3)), ',', num2str(x_best(4)), ',', num2str(x_best(5)), ',', num2str(x_best(6)), ',', num2str(x_best(7)), '] with an error of ', num2str(cost_best)])


