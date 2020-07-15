function [par_fit, cost_total] = fit_model1_Qmax_nonumves()

%Fit single-sensor model, fitting Qmax, not correcting num_ves.

Qmax_init = 6;

costfunc = @(Q_max)cost_model1_Q_max_nonumves(Q_max, 0);
options = optimset('Display', 'off', 'MaxIter', 100000, 'MaxFunEvals',10000);

tic;
[par_fit, cost_total] = fminsearch(costfunc, Qmax_init, options);
elaps_tim = toc




[cost_value, peak1s] = cost_model1_Q_max_nonumves(par_fit, 1);

save('./Sim_data/Log_files/fitresult_Q_max_notvescorrected.mat', 'par_fit', 'cost_total', 'peak1s')


save('./Sim_data/Log_files/best_fit_mod1_no_krep_no_numves.mat', 'par_fit', 'num_ves_factor')


filecont1 = sprintf('Optimization of model 1: One sensor. Fitting only Q_max. \n');

filecont2 = sprintf(     ['Q_max: ' num2str(par_fit) '\n' ...
                            'Cost: ' num2str(cost_total) '\n'   ...
                            'Elapsed time:' num2str(elaps_tim) ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

fullfile = './Sim_data/Log_files/Optimization_result_model1_Qmax_nonumves';

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
