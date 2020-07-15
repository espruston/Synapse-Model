function [cost_value, num_ves_factor, peak1s_corrected, peak2s_extr_corrected, pprs_extr] = cost_model2_Qmax_kon_s_corrected_numves(pars, save_data)

%%%This script determines the cost value of dual-sensor model.
%%%The num_ves is corrected for each choice of [Q_max, k_on, s] values to minimize the
%%%cost function

Q_max = pars(1);
k_on = pars(2);
s_val = pars(3);

CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;

disp('Running next simulation')
disp(['Q_max = ' num2str(Q_max)])
disp(['k_on = ' num2str(k_on)])
disp(['s_val = ' num2str(s_val)])
disp('Number of vesicles will be corrected later')

tic

%%%The Log file is written in this script

stoch_on_off = 0;
rand_ves_on_off = 1;
CalC_on_off = 1;
model_type = 24;


if Q_max <= 0
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_corrected = NaN;
else
    
    
    SS_coop = 2;
    k_D_second = 1.5e-6; %%k_D value is from Brandt et al. 2012
    
    
    par_free = [Q_max SS_coop k_D_second k_on s_val];
%     model_type = 24;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
    
[~, ~, ~, peak1s, peak2s_extr, pprs_extr, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc); 
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
                             'SS_ccop: ' num2str(SS_coop) '\n' ...
                             'k_on: ' num2str(k_on) '\n' ...
                             's_val: ' num2str(s_val) '\n' ...
                             'Num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                             'Cost: ' num2str(cost_value) '\n' ...
                            '------------------- \n', ...
                            'PVR: ' num2str(pVr1s) '\n' ...
                            'peak1s' num2str(peak1s_corrected) '\n', ...
                            'peak2s_extr' num2str(peak2s_extr) '\n', ...
                            'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

    fullfile = ['./Sim_data/Log_files/Log_Optimization_model2_Qmax_kon_s_numvesfactor'];

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');




toc