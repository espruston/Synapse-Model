function [cost_value, num_ves_factor, peak1s_corrected, peak2s_extr_corrected, pprs_extr] = cost_model1_Qmax_krep_corrected_numves(pars, save_data)

%%%This script determines the cost value of single-sensor model.
%%%The num_ves is corrected for each choice of Q_max and k_rep value to minimize the
%%%cost function


Q_max = pars(1);
k_rep = pars(2);

CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;


disp('Running next simulation')
disp(['Qmax = ' num2str(Q_max)])
disp('Number of vesicles will be corrected later')

tic;

%%%The Log file is written in this script

stoch_on_off = 0;
rand_ves_on_off = 1;

CalC_on_off = 1;

if Q_max <= 0
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_corrected = NaN;
else
    
    par_free = [Q_max k_rep];
    model_type = 15;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
    
[~, ~, ~, peak1s, peak2s_extr, pprs_extr, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc); 
% toc                                               
% 
% peak1s = peak1s_init*1e9;
% peak2s = peak2s_init*1e9;


[cost_value, num_ves_factor] = determine_cost_bothpeaks_numves_factor(peak1s, peak2s_extr)

peak1s_corrected = peak1s*num_ves_factor;
peak2s_extr_corrected = peak2s_extr*num_ves_factor;

end

cost_value

%%%Write log file for each iteration


    filecont = sprintf(     ['Q_max: ' num2str(Q_max) '\n' ...
                             'k_rep: ' num2str(k_rep) '\n' ...
                             'Num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                             'Cost: ' num2str(cost_value) '\n' ...
                            '------------------- \n', ...
                            'PVR: ' num2str(pVr1s) '\n' ...
                            'peak1s' num2str(peak1s) '\n', ...
                            'peak2s_extr' num2str(peak2s_extr) '\n', ...
                            'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

    fullfile = ['./Sim_data/Log_files/Log_Optimization_model1_Qmax_krep_numvesfactor'];

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');


toc
