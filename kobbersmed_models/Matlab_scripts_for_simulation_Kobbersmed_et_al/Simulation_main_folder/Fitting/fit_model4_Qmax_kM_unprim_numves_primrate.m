function [par_fit, cost_total] = fit_model4_Qmax_kM_unprim_numves_primrate(prim_coop)

%Fit unpriming model, fitting Qmax, kM_prim, unprim_rate, num_ves.

Qmax_init = 13.3025;
% Ca_prim_type = 5;
prim_kM_init = 4.8853e-08;
prim_rate_const_init = 181.82;
unprim_rate_const_init = 320.9958;

pars_init = [Qmax_init prim_kM_init prim_rate_const_init unprim_rate_const_init];

costfunc = @(pars)cost_model4_Qmax_kM_unprim_numves_primrate(pars, prim_coop, 0);
options = optimset('Display', 'off', 'MaxIter', 1000000, 'MaxFunEvals',100000);

tic;
[par_fit, cost_total] = fminsearch(costfunc, pars_init, options);
elaps_tim = toc


[cost_value, num_ves_factor, peak1s_corrected, peak2s_extr, pprs_extr] = cost_model4_Qmax_kM_unprim_numves_primrate(par_fit, prim_coop, 1);

save('./Sim_data/Log_files/fitresult_Q_max_numvescorrected.mat', 'par_fit', 'cost_value', 'num_ves_factor', 'peak1s_corrected', 'peak2s_extr')

save('./Sim_data/Log_files/best_fit_mod4_primrate.mat', 'par_fit', 'num_ves_factor')


filecont1 = sprintf('Optimization of model 4: One sensor, unpriming. Fitting Q_max, kM, unprim_rate, prim_rate and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['***Optimized parameters***' '\n' ...
                          'Q_max: ' num2str(par_fit(1)) '\n' ...
                          'prim_coop: ' num2str(prim_coop) '\n' ...
                          'k_M_prim: ' num2str(par_fit(2)) '\n' ...
                          'prim_const: ' num2str(par_fit(3)) '\n' ...
                          'unprim_const: ' num2str(par_fit(4)) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                          '***RESULTS***' '\n' ... 
                          'peak1s' num2str(peak1s_corrected) '\n', ...
                          'peak2s_extr' num2str(peak2s_extr) '\n', ...
                          'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                          'Cost: ' num2str(cost_value) '\n'   ...
                          'Elapsed time:' num2str(elaps_tim) ...
                          '\n' '------------------- \n' '\n' ...
                          ]);

fullfile = ['./Sim_data/Log_files/Optimization_result_model4_Qmax_kM_unprimrate_numvescorrected_primrate_primcoop' num2str(prim_coop)];

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
