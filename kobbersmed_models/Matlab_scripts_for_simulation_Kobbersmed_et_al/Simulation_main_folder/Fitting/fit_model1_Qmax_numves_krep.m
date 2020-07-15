function [par_fit, cost_total] = fit_model1_Qmax_numves_krep()

%Fit single-sensor model, fitting Qmax, krep, and num_ves.

Qmax_init = 6;
k_rep_init = 180;

pars_init = [Qmax_init k_rep_init];

costfunc = @(pars)cost_model1_Qmax_krep_corrected_numves(pars, 0);
options = optimset('Display', 'off', 'MaxIter', 100000, 'MaxFunEvals',10000);

tic;
[par_fit, cost_total] = fminsearch(costfunc, pars_init, options);
elaps_tim = toc



[cost_value, num_ves_factor, peak1s_corrected, peak2s_extr_corrected, pprs_extr] = cost_model1_Qmax_krep_corrected_numves(par_fit, 1);

save('./Sim_data/Log_files/fitresult_Q_max_numvescorrected.mat', 'par_fit', 'cost_total', 'num_ves_factor', 'peak1s_corrected')

save('./Sim_data/Log_files/best_fit_mod1_krep.mat', 'par_fit', 'num_ves_factor')


filecont1 = sprintf('Optimization of model 1: One sensor. Fitting Q_max, k_rep and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['Q_max: ' num2str(par_fit(1)) '\n' ...
                          'k_rep: ' num2str(par_fit(2)) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                          'Cost: ' num2str(cost_total) '\n'   ...
                          'Elapsed time:' num2str(elaps_tim) ...
                          '\n' '------------------- \n' '\n' ...
                          ]);

fullfile = './Sim_data/Log_files/Optimization_result_model1_Qmax_numvescorrected_krep';

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
