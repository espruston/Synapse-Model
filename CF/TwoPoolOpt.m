x0 = [0.57676,0.004604,0.00118,.75,1]; %[p_release, k_docking, k_undocking, size_rel,C_3]
lb = [.4, 0.0001, 0.0001,0.01,1];
ub = [1, 0.01, 0.01,1,10];

opts = optimoptions('surrogateopt','CheckpointFile','checkfile2pool.mat','PlotFcn','surrogateoptplot','UseParallel',true,'MaxTime',6*3600);
[x,err,exitflag,output] = surrogateopt(@TwoPoolFunc,lb,ub,opts);
%[x,err,exitflag,output] = surrogateopt('checkfile2pool.mat',opts);

disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)),'] with an error of ', num2str(err)])