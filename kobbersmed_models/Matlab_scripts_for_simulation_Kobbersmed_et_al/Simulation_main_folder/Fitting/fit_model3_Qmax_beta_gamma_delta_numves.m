function [par_fit, cost_total] = fit_model3_Qmax_beta_gamma_delta_numves()

%Fit site activation model, fitting Qmax, beta, gamma, delta, k_rep, num_ves.

Qmax_init = 13.3025;
% Ca_prim_type = 5;
beta_init = 0.12;
gamma_init = 150;
delta_init = 10;
k_rep_init = 180;

pars_init = [Qmax_init beta_init gamma_init delta_init k_rep_init];

costfunc = @(pars)cost_model3_Qmax_beta_gamma_delta_numves(pars, 0);
options = optimset('Display', 'off')%, 'MaxIter', 1000000, 'MaxFunEvals',100000);

tic;
[par_fit, cost_total] = fminsearch(costfunc, pars_init, options);
elaps_tim = toc


[cost_value, num_ves_factor, peak1s_corrected, peak2s_extr, pprs_extr] = cost_model3_Qmax_beta_gamma_delta_numves(par_fit, 1);

save('./Sim_data/Log_files/fitresult_mod3_krep_numves.mat.mat', 'par_fit', 'cost_value', 'num_ves_factor', 'peak1s_corrected', 'peak2s_extr')

save('./Sim_data/Log_files/best_fit_mod3_krep_numves.mat', 'par_fit', 'num_ves_factor')


filecont1 = sprintf('Optimization of model 3: One sensor, activation. Fitting Q_max, beta, gamma, delta, k_rep and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['***Optimized parameters***' '\n' ...
                          'Q_max: ' num2str(par_fit(1)) '\n' ...
                          'beta: ' num2str(par_fit(2)) '\n' ...
                          'gamma: ' num2str(par_fit(3)) '\n' ...
                          'delta: ' num2str(par_fit(4)) '\n' ...
                          'k_rep: ' num2str(par_fit(5)) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                          '***RESULTS***' '\n' ... 
                          'peak1s' num2str(peak1s_corrected) '\n', ...
                          'peak2s_extr' num2str(peak2s_extr) '\n', ...
                          'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                          'Cost: ' num2str(cost_value) '\n'   ...
                          'Elapsed time:' num2str(elaps_tim) ...
                          '\n' '------------------- \n' '\n' ...
                          ]);

fullfile = './Sim_data/Log_files/Optimization_result_model3_Qmax_beta_gamma_delta_numves';

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
