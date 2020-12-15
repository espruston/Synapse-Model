x0 = [.27; .071; 10; 0.015; .00012; 0.000023]; %initial guess p_mature = x(1); p_immature = x(2); C_3 = x(3); k_docking = x(4); k_maturation = x(5); k_dematuration = x(6);
A = [-1 1 0 0 0 0; 0 0 0 0 -1 1]; %p_immature - p_mature =< 0, k_dematuration - k_maturation =< 0  
b = [0; 0];
Aeq = [];
beq = [];
lb = [.05, .01, 0, .0001, .00001, 0];
ub = [1, .8, 100, 1, 1, 1];
nonlcon = [];
options = optimoptions('patternsearch','UseParallel',true,'MaxTime',3600);

[x,err] = patternsearch(@MaturationCDRFunc,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
disp(['Best fit was [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), '] with an error of ', num2str(err)])