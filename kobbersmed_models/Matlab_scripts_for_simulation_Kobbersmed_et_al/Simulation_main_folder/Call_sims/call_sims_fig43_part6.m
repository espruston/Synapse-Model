function call_sims_fig43_part6()
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
% CalC_on_off = 1;

% stoch_on_off = 1;
rand_ves_on_off = 1;

par_free(1) = 13.3025;
par_free(2) = 5;
par_free(3) = 4.8853e-08;
par_free(4) = 181.82;
par_free(5) = 320.9958;
par_free(6) = round(180*0.9376);
model_type = 42;

save_data = 2;

CaExtracellular = [8 8.5 9 9.5 10];

num_calc = length(CaExtracellular);

[par_init, savefilename] = parameter_choices(par_free, model_type);


testing_the_system_many_calcium(999, rand_ves_on_off, 1, par_init, save_data, savefilename, CaExtracellular);



for k = 1:num_calc
    CaExt = CaExtracellular(k);
    testing_the_system_morereps(1, rand_ves_on_off, par_init, save_data, savefilename, CaExt);
end

toc
