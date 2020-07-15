function [cost_value, peak1s] = cost_model1_Q_max_nonumves(Q_max, save_data)

%%%This script determines the cost value of single-sensor model.
%%%Q_max is optimized (num_ves is not). 


disp('Running next simulation')
disp(['Qmax = ' num2str(Q_max)])

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
    
    par_free = Q_max;
    model_type = 11;
    [par_init, savefilename] = parameter_choices(par_free, model_type);
    
tic; 
[~, ~, ~, peak1s, peak2s_extr, pprs_extr, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename); 
% toc                                               
% 
% peak1s = peak1s_init*1e9;
% peak2s = peak2s_init*1e9;


[cost_value] = determine_cost_peak1(peak1s);


end

cost_value

%%%Write log file for each iteration


    filecont = sprintf(     ['Q_max: ' num2str(Q_max) '\n' ...
                             'Cost: ' num2str(cost_value) '\n' ...
                            '------------------- \n', ...
                            'PVR: ' num2str(pVr1s) '\n' ...
                            'peak1s' num2str(peak1s) '\n', ...
                            'peak2s_extr' num2str(peak2s_extr) '\n', ...
                            'PPRs_extr: ' num2str(pprs_extr) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

    fullfile = ['./Sim_data/Log_files/Log_Optimization_Fig1_Qmax'];

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');



