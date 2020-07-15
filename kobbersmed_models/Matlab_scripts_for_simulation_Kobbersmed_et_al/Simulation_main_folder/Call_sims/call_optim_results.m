function call_optim_results(model_type, par_free, Ca_choice, stim_freq, num_stim)
%This file can be modified whenever the parameters need to be changed. par_init is used for optimisations.
%Model 11: par_free = Q_max
%Model 12: par_free = [Q_max kM]
%Model 13: par_free = [Q_max], k_rep = 0
%Model 14: par_free = [Q_max, num_ves], k_rep = 181.82
%Model 15: par_free = [Q_max, k_rep]
%Model 16: par_free = [Q_max, k_rep, num_ves]
%Model 21: par_free = [Q_max, SS_coop, s]. SS_PM = 0.
%Model 22: par_free = [Q_max, SS_coop, s]. SS_PM = 1.
%Model 23: par_free = [Q_max, SS_coop, k_d_second, k_4, s]. SS_PM = 0
%Model 24: par_free = [Q_max, SS_coop, k_d_second, k_4, s]. SS_PM = 1
%Model 231: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep]. SS_PM = 0
%Model 241: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep]. SS_PM = 1
%Model 242: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep, num_ves]. SS_PM = 1
%Model 25: par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 0
%Model 26: par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 1
%Model 27: par_free = [Q_max, SS_coop, s]. SS_PM = 0, k_rep = 0;
%Model 28: par_free = [Q_max, SS_coop, s]. SS_PM = 1, k_rep = 0;
%Model 31: par_free = [Q_max, act_model_type, beta, gamma, delta]
%Model 32: par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep]
%Model 33: par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep, num_ves]
%Model 41: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
%Model 42: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
%Model 43: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, k_M_Rest, CaMax_rest]
%Model 44: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0]
%Model 45: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, CaMax_rest]
%Model 51: par_free = [Q_max, SS_coop, prim_rate, unprim_rate, u]; %Regehr values, SS_PM = 0.

CalC_on_off = 1;
stoch_on_off = 1;
rand_ves_on_off = 1;

% save_data = 1; %1: a lot, 2: less, 66: only individual peaks and EPSCs
save_calc_loc = 1; %1 for usual (./Sim_data/CalC_files), 2 for alternative (./Sim_data/new_calcium)

num_iterations = 10; %call 10 sets of 100 simulations
% num_iterations = 5 %call 5 sets of 100 simulations
% disp('Simulating only 500 reps')


if Ca_choice == 1
    CaExtracellular = [0.75 1 1.25 1.5 1.75 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]; %All concentrations
    disp('Simulating a lot of calcium concentrations')
elseif Ca_choice == 2
    CaExtracellular = [ 1 1.25 1.75 2 2.5 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5]; %All concentrations except experimental
    disp('Simulating a lot of calcium - except experimental concentrations')
elseif Ca_choice == 3
    CaExtracellular = [0.75 1.5 3 6 10]; %Experimental concentrations
    disp('Only simulating experimental concentrations')
elseif Ca_choice == 4
    CaExtracellular = [8 8.5 9 9.5 10]; %Some arbitrary concentrations, changed depending on purpose
    disp('Only simulating some concentrations')
end

num_calc = length(CaExtracellular);


% %Deterministic
% if stoch_on_off == 2 || stoch_on_off == 0
%     save_data = 1;
% [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
% testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc);
% CalC_on_off = 0;
% end



%Stochastic

    for k = 1:num_calc
%         tic;
        run_more_reps(num_iterations, par_free, model_type, CalC_on_off, stoch_on_off, rand_ves_on_off, CaExtracellular(k), save_calc_loc, stim_freq, num_stim)
%         toc;
    end





