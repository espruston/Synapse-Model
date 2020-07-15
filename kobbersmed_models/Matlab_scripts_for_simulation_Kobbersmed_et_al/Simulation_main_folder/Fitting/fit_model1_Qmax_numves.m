function [par_fit, cost_total] = fit_model1_Qmax_numves()

%Fit single-sensor model, fitting Qmax and num_ves.

Qmax_init = 6;

costfunc = @(Q_max)cost_model1_Qmax_corrected_numves(Q_max, 0);
options = optimset('Display', 'off', 'MaxIter', 100000, 'MaxFunEvals',10000);

tic;
[par_fit, cost_total] = fminsearch(costfunc, Qmax_init, options);
elaps_tim = toc



[cost_value, num_ves_factor, peak1s, peak1s_corrected] = cost_model1_Qmax_corrected_numves(par_fit, 1);

save('./Sim_data/Log_files/fitresult_Q_max_numvescorrected.mat', 'par_fit', 'cost_total', 'num_ves_factor', 'peak1s_corrected')

save('./Sim_data/Log_files/best_fit_mod1_no_krep.mat', 'par_fit', 'num_ves_factor')



filecont1 = sprintf('Optimization of model 1: One sensor. Fitting Q_max and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['Q_max: ' num2str(par_fit) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                            'Cost: ' num2str(cost_total) '\n'   ...
                            'Elapsed time:' num2str(elaps_tim) ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

fullfile = './Sim_data/Log_files/Optimization_result_model1_Qmax_numvescorrected';

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
