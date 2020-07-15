function call_calcium_simulation()
%This file can be modified whenever the parameters need to be changed. par_init is used for optimisations.
%Model 11: par_free = Q_max
%Model 12: par_free = [Q_max kM]
%Model 13: par_free = [Q_max], k_rep = 0
%Model 14: par_free = [Q_max, num_ves], k_rep = 181.82
%Model 21: par_free = [Q_max, SS_coop, s]. SS_PM = 0.
%Model 22: par_free = [Q_max, SS_coop, s]. SS_PM = 1.
%Model 23: par_free = [Q_max, SS_coop, k_4, k_min4, s]. SS_PM = 0
%Model 24: par_free = [Q_max, SS_coop, k_4, k_min4, s]. SS_PM = 1
%Model 31: par_free = [Q_max, act_model_type, kMact, gamma, delta]
%Model 41: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
%Model 42: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
%Model 51: par_free = [Q_max, SS_coop, prim_rate, unprim_rate, u]; %Regehr values, SS_PM = 0.
%Model 51: par_free = [Q_max, SS_coop, k_4, k_min4, u]. SS_PM = 0

% 




%%%CALCIUM


CaExtracellular = [0.75 1.5 3 6 10];

rand_ves_on_off = 0;
save_data = 0;
save_calc_loc = 2; 
pVr2_hack = 0;
stoch_on_off = 999;
CalC_on_off = 999;

Q_maxes = [8.4201 4.5169 12.5872 13.7718];

for k = 1:length(Q_maxes)
    
    Q_max = Q_maxes(k);
    


    disp(['CalC simulation.'])

    par_free = Q_max;
    model_type = 11;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);
    
    tic;
%     testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
    disp('CalC simulation time:')
    toc
    
end
