function [par_fit, cost_total] = fit_model2_Qmax_kon_s_numves_krep_perimeter(SS_coop)

%Fit dual-sensor model (SVs at single distance), fitting Qmax, k_on, s, num_ves, k_rep.

Qmax_init = 6;
k_on_init = 1e2; %Will be multiplied with 1e4
s_val_init = 200;
k_rep_init = 180;

pars_init = [Qmax_init k_on_init s_val_init k_rep_init];

costfunc = @(pars)cost_model2_Qmax_kon_s_corrected_numves_krep_perimeter(pars, 0, SS_coop);
options = optimset('Display', 'off');

tic;
[par_fit, cost_total] = fminsearch(costfunc, pars_init, options);
elaps_tim = toc



[cost_value, num_ves_factor, peak1s_corrected, peak2s_extr_corrected, pprs_extr] = cost_model2_Qmax_kon_s_corrected_numves_krep_perimeter(par_fit, 1, SS_coop);

save(['./Sim_data/Log_files/best_fit_mod2_SScoop' num2str(SS_coop) '_krep_Q_max_numvescorrected_perimeter.mat'], 'par_fit', 'num_ves_factor', 'peak1s_corrected', 'peak2s_extr_corrected', 'cost_total')



filecont1 = sprintf('Optimization of model 2: Two sensors. Fitting Q_max, k_on, s, k_rep and a factor on the number of vesicles. \n');

filecont2 = sprintf(     ['***Optimized parameters***' '\n' ...
                          'Q_max: ' num2str(par_fit(1)) '\n' ...
                          'k_on: ' num2str(par_fit(2)) '\n' ...
                          's_val: ' num2str(par_fit(3)) '\n' ...
                          'k_rep: ' num2str(par_fit(4)) '\n' ...
                          'num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                          '***RESULTS***' '\n' ... 
                          'peak1s' num2str(peak1s_corrected) '\n', ...
                          'peak2s_extr' num2str(peak2s_extr_corrected) '\n', ...
                          'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                          'Cost: ' num2str(cost_total) '\n'   ...
                          'Elapsed time:' num2str(elaps_tim) ...
                          '\n' '------------------- \n' '\n' ...
                          ]);

fullfile = ['./Sim_data/Log_files/Optimization_result_model2_SScoop' num2str(SS_coop) '_Qmax_kon_s_numvescorrected_krep_perimeter'];

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont1);
fprintf(newParameterfile,'%s\n',filecont2);
fclose(newParameterfile);

fclose('all');
