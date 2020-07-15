function [cost_value, num_ves_factor, peak1s_corrected, peak2s_extr_corrected, pprs_extr] = cost_model2_Qmax_kon_s_corrected_numves_krep(pars, save_data, SS_coop, k_D_second)

%%%This script determines the cost value of dual-sensor model.
%%%The num_ves is corrected for each choice of [Q_max, k_on, s, k_rep] values to minimize the
%%%cost function

Q_max = pars(1);
k_on = pars(2)*1e4;
s_val = pars(3);
k_rep = pars(4);

CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;

disp('Running next simulation')
disp(['Q_max = ' num2str(Q_max)])
disp(['k_on = ' num2str(k_on)])
disp(['s_val = ' num2str(s_val)])
disp(['k_rep = ' num2str(k_rep)])
disp('Number of vesicles will be corrected later')

%%%The Log file is written in this script

stoch_on_off = 0;
rand_ves_on_off = 1;
CalC_on_off = 1;
model_type = 241;
pVr2_hack = 0;

% if fig52_on_off == 1
%     CalC_on_off = 0;
% end


if Q_max <= 0
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_corrected = NaN;
else
    
    
%     SS_coop = 2;
%     k_D_second = 1.5e-6; %%k_D value is from Brandt et al. 2012

    
    
    par_free = [Q_max SS_coop k_D_second k_on s_val k_rep];
%     model_type = 24;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);
    
[~, ~, ~, peak1s, peak2s_extr, pprs_extr, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack); 
% 
% peak1s = peak1s_init*1e9;
% peak2s = peak2s_init*1e9;


[cost_value, num_ves_factor] = determine_cost_bothpeaks_numves_factor(peak1s, peak2s_extr);

peak1s_corrected = peak1s*num_ves_factor;
peak2s_extr_corrected = peak2s_extr*num_ves_factor;

end

num_ves_factor
cost_value

%%%Write log file for each iteration


    filecont = sprintf(     ['Q_max: ' num2str(Q_max) '\n' ...
                             'SS_coop: ' num2str(SS_coop) '\n' ...
                             'k_D_second: ' num2str(k_D_second) '\n' ...                             'k_on: ' num2str(k_on) '\n' ...
                             's_val: ' num2str(s_val) '\n' ...
                             'k_rep: ' num2str(k_rep) '\n' ...
                             'Num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                             'Cost: ' num2str(cost_value) '\n' ...
                            '------------------- \n', ...
                            'PVR: ' num2str(pVr1s) '\n' ...
                            'peak1s' num2str(peak1s_corrected) '\n', ...
                            'peak2s_extr' num2str(peak2s_extr_corrected) '\n', ...
                            'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);


% if fig52_on_off == 1    
%         fullfile = './Sim_data/Log_files/Log_fig52_simulation';
% 
%         file_pars_string = [];
%         for kk = 1:length(pars)
%             file_pars_string = [file_pars_string '_' num2str(pars(kk))];
%         end
% 
%         save(['./Sim_data/new_results/Fig52_results_pars' file_pars_string '.mat'], 'peak1s_corrected', 'peak2s_extr_corrected', 'pprs_extr', 'pVr1s', 'pars', 'k_on', 's_val', 'k_rep', 'SS_coop', 'num_ves_factor', 'Q_max', 'cost_value',  '-v7.3')
% 
% else
    fullfile = ['./Sim_data/Log_files/Log_Optimization_model2_SScoop' num2str(SS_coop) '_Qmax_kon_s_numvesfactor_krep_kD_' num2str(k_D_second)];
% end

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');


