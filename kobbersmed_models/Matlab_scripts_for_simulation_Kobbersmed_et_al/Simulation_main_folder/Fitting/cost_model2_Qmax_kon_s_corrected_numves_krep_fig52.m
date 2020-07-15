function [cost_value_peak1, cost_value_bothpeaks, num_ves_factor_peak1, num_ves_factor_bothpeaks, peak1s_correctedtopeak1, peak1s_correctedtobothpeaks, peak2s_extr_correctedtopeak1, peak2s_extr_correctedtobothpeaks, pprs_extr] = cost_model2_Qmax_kon_s_corrected_numves_krep_fig52(pars, save_data)

%%%This script determines the cost value of dual-sensor model.
%%%The num_ves is corrected for each choice of [Q_max, k_on, s, k_rep] values to minimize the
%%%cost function

%%%NOTE: This file is used when generating files for plot 52 (3d plots in fig. 6)

Q_max = pars(1);
k_D_second = pars(2);
k_on = pars(3);
s_val = pars(4);
k_rep = pars(5);

CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;

disp('Running next simulation')
disp(['Q_max = ' num2str(Q_max)])
disp(['k_D = ' num2str(k_D_second)])
disp(['k_on = ' num2str(k_on)])
disp(['s_val = ' num2str(s_val)])
disp(['k_rep = ' num2str(k_rep)])
disp('Number of vesicles will be corrected later')

%%%The Log file is written in this script

stoch_on_off = 0;
rand_ves_on_off = 1;
% CalC_on_off = 1;
model_type = 241;


% if fig52_on_off == 1
    CalC_on_off = 0;
% end


if Q_max <= 0
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_correctedtopeak1 = NaN;
else
    
    
    SS_coop = 2;
%     k_D_second = 1.5e-6; %%k_D value is from Brandt et al. 2012

    
    
    par_free = [Q_max SS_coop k_D_second k_on s_val k_rep];
%     model_type = 24;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);
    
    pVr2_hack = 0;
    
[~, ~, ~, peak1s, peak2s_extr, pprs_extr, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack); 
% 
% peak1s = peak1s_init*1e9;
% peak2s = peak2s_init*1e9;


[cost_value_bothpeaks, num_ves_factor_bothpeaks] = determine_cost_bothpeaks_numves_factor(peak1s, peak2s_extr);
[cost_value_peak1, num_ves_factor_peak1] = determine_cost_peak1_numves_factor(peak1s);


peak1s_correctedtopeak1 = peak1s*num_ves_factor_peak1;
peak1s_correctedtobothpeaks = peak1s*num_ves_factor_bothpeaks;
peak2s_extr_correctedtopeak1 = peak2s_extr*num_ves_factor_peak1;
peak2s_extr_correctedtobothpeaks = peak2s_extr*num_ves_factor_bothpeaks;

end

num_ves_factor_peak1
cost_value_peak1

num_ves_factor_bothpeaks
cost_value_bothpeaks

%%%Write log file for each iteration


    filecont = sprintf(     ['Q_max: ' num2str(Q_max) '\n' ...
                             'SS_coop: ' num2str(SS_coop) '\n' ...
                             'k_on: ' num2str(k_on) '\n' ...
                             's_val: ' num2str(s_val) '\n' ...
                             'k_rep: ' num2str(k_rep) '\n' ...
                             'Num_ves_factor peak1: ' num2str(num_ves_factor_peak1) '\n' ...
                             'Num_ves_factor both peaks: ' num2str(num_ves_factor_bothpeaks) '\n' ...
                             'Cost, peak1: ' num2str(cost_value_peak1) '\n' ...
                             'Cost, both peak: ' num2str(cost_value_bothpeaks) '\n' ...
                            '------------------- \n', ...
                            'PVR: ' num2str(pVr1s) '\n' ...
                            'peak1s corrected to peak1' num2str(peak1s_correctedtopeak1) '\n', ...
                            'peak1s corrected to both peaks' num2str(peak1s_correctedtobothpeaks) '\n', ...
                            'peak2s_extr corrected to peak1' num2str(peak2s_extr_correctedtopeak1) '\n', ...
                            'peak2s_extr corrected to both peaks' num2str(peak2s_extr_correctedtobothpeaks) '\n', ...
                            'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);


% if fig52_on_off == 1    
        fullfile = './Sim_data/Log_files/Log_fig52_simulation';

        file_pars_string = [];
        for kk = 1:length(pars)
            file_pars_string = [file_pars_string '_' num2str(pars(kk))];
        end

        save(['./Sim_data/new_results/Fig52_results_pars' file_pars_string '.mat'], 'peak1s_correctedtopeak1', 'peak1s_correctedtobothpeaks', 'peak2s_extr_correctedtopeak1', 'peak2s_extr_correctedtobothpeaks', 'pprs_extr', 'pVr1s', 'pars', 'k_on', 's_val', 'k_rep', 'SS_coop', 'num_ves_factor_peak1', 'num_ves_factor_bothpeaks', 'Q_max', 'cost_value_peak1', 'cost_value_bothpeaks',  '-v7.3')

% else
%     fullfile = './Sim_data/Log_files/Log_Optimization_model2_Qmax_kon_s_numvesfactor_krep';
% end

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');


