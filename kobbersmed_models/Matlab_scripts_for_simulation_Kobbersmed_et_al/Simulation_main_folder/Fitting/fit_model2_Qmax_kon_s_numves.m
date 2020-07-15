function [par_fit, cost_total] = fit_model2_Qmax_kon_s_numves()

%Fit dual-sensor model, fitting Qmax, k_on, s, num_ves.

Qmax_init = 6;
k_on_init = 1e6;
s_val_init = 200;

pars_init = [Qmax_init k_on_init s_val_init];

costfunc = @(pars)cost_model2_Qmax_kon_s_corrected_numves(pars, 0);
options = optimset('Display', 'off', 'MaxIter', 1000000, 'MaxFunEvals',100000);

tic;
[par_fit, cost_total] = fminsearch(costfunc, pars_init, options);
elaps_tim = toc



[cost_value, num_ves_factor, peak1s_corrected, peak2s_extr, pprs_extr] = cost_model2_Qmax_kon_s_corrected_numves(par_fit, 1);

save('./Sim_data/Log_files/fitresult_Q_max_numvescorrected.mat', 'par_fit', 'cost_total', 'num_ves_factor', 'peak1s_corrected', 'peak2s_extr')

save('./Sim_data/Log_files/best_fit_mod2_no_krep.mat', 'par_fit', 'num_ves_factor')


filecont1 = sprintf('Optimization of model 2: Two sensors. Fitting Q_max, k_on, s and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['***Optimized parameters***' '\n' ...
                          'Q_max: ' num2str(par_fit(1)) '\n' ...
                          'k_on: ' num2str(par_fit(2)) '\n' ...
                          's_val: ' num2str(par_fit(3)) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                          '***RESULTS***' '\n' ... 
                          'peak1s' num2str(peak1s_corrected) '\n', ...
                          'peak2s_extr' num2str(peak2s_extr) '\n', ...
                          'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                          'Cost: ' num2str(cost_total) '\n'   ...
                          'Elapsed time:' num2str(elaps_tim) ...
                          '\n' '------------------- \n' '\n' ...
                          ]);

fullfile = './Sim_data/Log_files/Optimization_result_model2_Qmax_kon_s_numvescorrected';

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
