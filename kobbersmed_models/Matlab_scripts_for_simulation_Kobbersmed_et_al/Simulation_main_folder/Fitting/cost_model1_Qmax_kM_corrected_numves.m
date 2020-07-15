function [cost_value, num_ves_factor, peak1s, peak1s_corrected] = cost_model1_Qmax_kM_corrected_numves(Q_max_and_kM, save_data)

%%%This script determines the cost value of single-sensor model.
%%%The num_ves is corrected for each choice of Q_max and k_M value to minimize the
%%%cost function

disp('Running next simulation')
disp(['Qmax = ' num2str(Q_max_and_kM(1))])
disp(['kM = ' num2str(Q_max_and_kM(2))])
disp('Number of vesicles will be corrected later')

tic;
%%%The Log file is written in this script

stoch_on_off = 0;
k_rep = 181.8182;
act_model_type = 0;
k_M_act = 0;
rand_ves_on_off = 0;
gamma = 0;
delta = 0;
SS_coop = 0;
opt_s = 1;
opt_k_d_second = 0;
opt_k_4_second = 0;
num_ves = 180;
% save_data = 0;

Ca_prim_type = 0;
prim_kM = 0;
prim_rate_const = 0;
unprim_rate_const = 0;
SS_PM = 0;

CalC_on_off = 1;

if sum(Q_max_and_kM <= 0)
    cost_value = 1e6;
    num_ves_factor = NaN;
    peak1s = NaN;
    peak1s_corrected = NaN;
else


[~, ~, ~, peak1s, peak2s, pprs, ~, pVr1s, ~, ~, ~] = testing_the_system(stoch_on_off, act_model_type, Ca_prim_type, SS_PM, SS_coop, k_M_act, rand_ves_on_off, CalC_on_off, Q_max_and_kM, num_ves, k_rep, gamma, delta, prim_kM, prim_rate_const, unprim_rate_const, opt_s, opt_k_d_second, opt_k_4_second, save_data); 

% 
% peak1s = peak1s_init*1e9;
% peak2s = peak2s_init*1e9;


[cost_value, num_ves_factor] = determine_cost_numves_factor(peak1s);

peak1s_corrected = peak1s*num_ves_factor;

end

cost_value

%%%Write log file for each iteration


    filecont = sprintf(     ['Q_max: ' num2str(Q_max_and_kM(1)) '\n' ...
                             'k_M: ' num2str(Q_max_and_kM(2)) '\n' ...
                             'Num_ves_factor: ' num2str(num_ves_factor) '\n' ...
                             'Cost: ' num2str(cost_value) '\n' ...
                            '------------------- \n', ...
                            'PVR: ' num2str(pVr1s) '\n' ...
                            'peak1s' num2str(peak1s) '\n', ...
                            'peak2s' num2str(peak2s) '\n', ...
                            'PPR: ' num2str(pprs) '\n'   ...
                            '\n' '------------------- \n' '\n' ...
                            ]);

    fullfile = ['./Sim_data/Log_files/Log_Optimization_model1_Qmax_kM_numvesfactor'];

    newParameterfile = fopen(fullfile,'a');
    fprintf(newParameterfile,'%s\n',filecont);
    fclose(newParameterfile);

    fclose('all');


toc
