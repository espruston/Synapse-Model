function [cost_value_peak1, cost_value_bothpeaks, num_ves_factor_peak1, num_ves_factor_bothpeaks, peak1s_correctedtopeak1, peak1s_correctedtobothpeaks, peak2s_extr_correctedtopeak1, peak2s_extr_correctedtobothpeaks, pprs_extr] = cost_model4_fig54(pars, save_data)

%%%This script determines the cost value of unpriming model.
%%%The num_ves is corrected for each choice of [Q_max, k_M_prim, prim_rate, unprim_rate] values to minimize the
%%%cost function

%%NOTE: Only used when generating results for fig. 54 (3D plots in fig. 7)

Q_max = pars(1);
k_M_prim = pars(2);
prim_rate = pars(3);
unprim_rate = pars(4);


CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;

disp('Running next simulation')
disp(['Q_max = ' num2str(Q_max)])
disp(['k_M_prim = ' num2str(k_M_prim)])
disp(['prim_rate = ' num2str(prim_rate)])
disp(['unprim_rate = ' num2str(unprim_rate)])
disp('Number of vesicles will be corrected later')

%%%The Log file is written in this script

stoch_on_off = 0;
rand_ves_on_off = 1;
% CalC_on_off = 1;
model_type = 41;


% if fig52_on_off == 1
    CalC_on_off = 0;
% end


if Q_max <= 0
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_correctedtopeak1 = NaN;
else
    
    
    prim_type = 5;

    
    
    par_free = [Q_max prim_type k_M_prim prim_rate unprim_rate];
%     model_type = 24;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
    
[~, ~, ~, peak1s, peak2s_extr, pprs_extr, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc); 
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
                             'Prim_type: ' num2str(prim_type) '\n' ...
                             'k_M_prim: ' num2str(k_M_prim) '\n' ...
                             'prim_rate: ' num2str(prim_rate) '\n' ...
                             'unprim_rate: ' num2str(unprim_rate) '\n' ...
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
        fullfile = './Sim_data/Log_files/Log_mod4_fig54_simulation';

        file_pars_string = [];
        for kk = 1:length(pars)
            file_pars_string = [file_pars_string '_' num2str(pars(kk))];
        end

        save(['./Sim_data/new_results/Fig54_mod4_results_pars' file_pars_string '.mat'], 'peak1s_correctedtopeak1', 'peak1s_correctedtobothpeaks', 'peak2s_extr_correctedtopeak1', 'peak2s_extr_correctedtobothpeaks', 'pprs_extr', 'pVr1s', 'pars', 'unprim_rate', 'k_M_prim', 'prim_rate', 'prim_type', 'num_ves_factor_peak1', 'num_ves_factor_bothpeaks', 'Q_max', 'cost_value_peak1', 'cost_value_bothpeaks',  '-v7.3')

% else
%     fullfile = './Sim_data/Log_files/Log_Optimization_model2_Qmax_kon_s_numvesfactor_krep';
% end

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');


