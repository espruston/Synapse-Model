function call_single_simulation() %(model_type, par_free)
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
%Model 31: par_free = [Q_max, act_model_type, kMact, gamma, delta]
%Model 41: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
%Model 42: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
%Model 43: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, k_M_Rest, CaMax_rest]
%Model 44: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0]
%Model 45: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, CaMax_rest]
%Model 46: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
%Model 51: par_free = [Q_max, SS_coop, prim_rate, unprim_rate, u]; %Regehr values, SS_PM = 0.


% 
% % Model 4
% model_type = 46;
% par_free(1) = 13.7718;
% par_free(2) = 5;
% par_free(3) = 5.5207e-08;
% par_free(4) = 134.8542;
% par_free(5) = 236.8238;
% par_free(6) = 180;
% 
% stoch_on_off = 1;
% rand_ves_on_off = 1;
% CalC_on_off = 1;
% CaExtracellular = [0.75 1.5 3 6 10];
% save_data = 1;
% save_calc_loc = 1;
% pVr2_hack = 0;
% 
% [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);
% 
% % testing_the_system(0, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, 1, savefilename, save_calc_loc, pVr2_hack, 2, 200)
% 
% % testing_the_system(stoch_on_off, rand_ves_on_off, 0, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack, 100, 2)
% 
% 
% % for Ca = [3 6 10]
% %     testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, Ca, save_data, savefilename, save_calc_loc, pVr2_hack, 20, 2)
% % end
% 
% % ISIs = [100 250 500 1000 5 10 25 50];
% 
% ISIs = [250 500 1000 5 10 25 50];
% 
% ISIs = [10]
% 
% stim_freqs = 1000./ISIs;
% 
% num_stim = 2;
% 
% for stim_freq = stim_freqs
%     for Ca = CaExtracellular
%         testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, Ca, save_data, savefilename, save_calc_loc, pVr2_hack, stim_freq, num_stim)
%     end
% end
% stoch_on_off = 0





% Model 3
model_type = 33;

par_free = [12.587204133650637 5 0.092900455549552 1.947727198548635e+02 10.697511762243380 1.412098817171219e+02 round(1.051855468750000*180)];

stoch_on_off = 1;
rand_ves_on_off = 1;
CalC_on_off = 0;
CaExtracellular = [0.75 1.5 3 6 10];
save_data = 1;
save_calc_loc = 1;
pVr2_hack = 0;

num_stim = 2;


[par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);




ISIs = [250 500 1000 5 10 25 50 100];

ISIs = [100];

stim_freqs = 1000./ISIs;

stim_freqs = 4;

for stim_freq = stim_freqs
%     for Ca = [3 6 10]

            Ca = 3;
        
        testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, Ca, save_data, savefilename, save_calc_loc, pVr2_hack, stim_freq, num_stim)
%     end
end



