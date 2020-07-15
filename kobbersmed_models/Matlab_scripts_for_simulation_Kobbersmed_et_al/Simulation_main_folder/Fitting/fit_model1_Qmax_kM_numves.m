function [par_fit, cost_total] = fit_model1_Qmax_kM_numves()

%Fit single-sensor model, fitting Qmax, kM, and num_ves.



Qmax_init = 6;
kM_init = 2.679;

costfunc = @(Q_max_and_kM)cost_model1_Qmax_kM_corrected_numves(Q_max_and_kM, 0);
options = optimset('Display', 'off', 'MaxIter', 150000, 'MaxFunEvals',150000);

par_init = [Qmax_init kM_init];

tic;
[par_fit, cost_total] = fminsearch(costfunc, par_init, options);
elaps_tim = toc

Q_max_and_kM = par_fit;


[cost_value, num_ves_factor, peak1s, peak1s_corrected] = cost_model1_Qmax_kM_corrected_numves(Q_max_and_kM, 1)

save('./Sim_data/Log_files/fitresult_Qmax_kM_numvescorrected.mat', 'par_fit', 'cost_total', 'peak1s')

save('./Sim_data/Log_files/best_fit_mod1_km_no_krep.mat', 'par_fit', 'num_ves_factor')


filecont1 = sprintf('Optimization of model 1: One sensor. Fitting Q_max, kM and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['Q_max: ' num2str(par_fit(1)) '\n' ...
                          'kM: ' num2str(par_fit(2)) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                            'Cost: ' num2str(cost_total) '\n'   ...
                            'Elapsed time:' num2str(elaps_tim) ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

fullfile = './Sim_data/Log_files/Optimization_result_model1_Qmax_kM_numvescorrected';

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
