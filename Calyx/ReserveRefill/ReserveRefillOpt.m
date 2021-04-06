%Docking Increase
% x0 = [0.28747,0.0010252,3.9746e-05,7.2955,4.1053e-05,0,0]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3, C_7]
% lb = [.1, 0.0001, 0.00001, 1, 0.000001, 0, 0];
% ub = [.5, 0.008, 0.001, 100, 0.0001, 0, 0];
% A = [];
% b = [];
% fun = @Syt3DockingIncreaseFunc;

%Docking Increase 2 State
% x0 = [0.671166,0.599586,0.00569202,0.00572867,0.137335,0.0176688,27.407,0.0832238,30.8695,30]; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3, C_7]
% lb = [0, 0, 0.0001, 0.0001, .0001, .00001, 10, 0.00001, 0, 0];
% ub = [1, 1, 0.1, 0.1, 0.9, 0.9, 50, 0.1, 50, 50];
% A = [];
% b = [];
% fun = @DockingIncrease2StateFunc;

%Docking Increase 2 State DKO
x0 = [0.4,0.7,0.01,0.0001,0.0002,5e-05,100,1e-05] ; %[p_release, k_docking, k_undocking, reserve_size, k_refill, C_3, C_7]
lb = [0, 0, 0.0001, 0.0001, .0001, .00001, 10, 0.00001];
ub = [1, 1, 0.5, 0.1, 0.5, 0.1, 200, 0.5];
A = [];
b = [];
fun = @DockingIncrease2StateDKOFunc;

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

user_input_time = 0.25*3600;
count = 0; count_max = 2000;
x = x0;
cost_best = 100;
x_best = x0;
tic;
while (toc < user_input_time) && (count < count_max)
    
    x = lb + rand(size(x0)).*(ub-lb);
    
    cost = fun(x);
    
    if cost < cost_best
        cost_best = cost;
        x_best = x;
    end
    
    count = count + 1;
end


options = optimoptions('patternsearch','PlotFcn','psplotbestf','MaxTime',1.5*3600);
[x_best,cost_best] = patternsearch(fun,x_best,[],[],[],[],lb,ub,[],options);

disp(['Best fit was [', num2str(x_best(1)), ',', num2str(x_best(2)), ',', num2str(x_best(3)), ',', num2str(x_best(4)), ',', num2str(x_best(5)), ',', num2str(x_best(6)), ',', num2str(x_best(7)), ',', num2str(x_best(8)), '] with an error of ', num2str(cost_best)])

%disp(['Best fit was [', num2str(x_best(1)), ',', num2str(x_best(2)), ',', num2str(x_best(3)), ',', num2str(x_best(4)), ',', num2str(x_best(5)), ',', num2str(x_best(6)), ',', num2str(x_best(7)), '] with an cost of ', num2str(cost_best)])
