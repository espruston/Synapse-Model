function call_sims_figUP()
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
tic;
CalC_on_off = 1;

stoch_on_off = 2;
rand_ves_on_off = 1;

par_free(1) = 14.5;
par_free(2) = 3;
par_free(3) = 5e-8;
par_free(4) = 150;
par_free(5) = 700;
model_type = 41;

save_data = 2;

[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_many_calcium(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
toc

