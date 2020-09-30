x0 = [.67; .2; .03; 0; .0001; 0]; %initial guess p_mature = x(1); p_immature = x(2); k_docking = x(3); k_undocking = x(4); k_maturation = x(5); k_dematuration = x(6);
A = [-1 1 0 0 0 0; 0 0 -1 1 0 0; 0 0 0 0 -1 1]; %p_immature - p_mature =< 0, k_undocking - k_docking =< 0, k_dematuration - k_maturation =< 0  
b = [0; 0; 0];
Aeq = [];
beq = [];
lb = [.1, .01, .0001, 0, .000001, 0];
ub = [1, .8, 1, 1, 1, 1];
nonlcon = [];
options = optimoptions('patternsearch','UseParallel',true);

x = patternsearch(@MaturationCDRFunc,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);