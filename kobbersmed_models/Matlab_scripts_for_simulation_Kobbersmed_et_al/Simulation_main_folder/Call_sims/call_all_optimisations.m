function call_all_optimisations()


% fit_model2_Qmax_kon_s_numves_krep(2, 1.5e-6)

% for coop = 4:5
%     fit_model2_Qmax_kon_s_numves_krep(coop, 1e-6)
%     fit_model2_Qmax_kon_s_numves_krep(coop, 2e-6)
% end
% 
% for coop = 3:5
%     fit_model2_Qmax_kon_s_numves_krep(coop, 5e-7)
%     fit_model2_Qmax_kon_s_numves_krep(coop, 3e-6)
% end




model_type = 46;



stoch_on_off = 1;
rand_ves_on_off = 1;
CalC_on_off = 1;
iteration = 0;
save_data = 3;
CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;
pVr2_hack = 0;

par_free = [13.7718 5 5.5207e-08 134.8542 236.8238 round(1.0025*180)];
[par_init, savefilename] = parameter_choices(par_free, model_type, iteration, rand_ves_on_off);


ISIs = [5 10 25 50 100 250 500 1000];

stim_freqs = 1000./ISIs;

num_stim = 2;

for stim_freq = stim_freqs
    call_optim_results(model_type, par_free, 3, stim_freq, num_stim)
%     call_optim_results(model_type, par_free, 2, stim_freq, num_stim)
end


