function call_all_simulations()


model_type = 46;



stoch_on_off = 1;
rand_ves_on_off = 1;
CalC_on_off = 1;
iteration = 0;
save_data = 3;
CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 2;
pVr2_hack = 0;

par_free = [13.7718 5 5.5207e-08 134.8542 236.8238 round(1.0025*180)];
[par_init, savefilename] = parameter_choices(par_free, model_type, iteration, rand_ves_on_off);
testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack)

par_free = [13.4941 2 7.612e-09 106.5916 5207.7013 round(1.1379*180)];
[par_init, savefilename] = parameter_choices(par_free, model_type, iteration, rand_ves_on_off);
testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack)



Ca_choice = 1;
% 
%%Model 4: Only one unpriming reaction
model_type = 46;

par_free = [13.7718 5 5.5207e-08 134.8542 236.8238 round(1.0025*180)];

call_optim_results(model_type, par_free, 4);


par_free = [13.4941 2 7.612e-09 106.5916 5207.7013 round(1.1379*180)];

call_optim_results(model_type, par_free, Ca_choice);

%%%%%%






% 
% Ca_choice = 1;
% % 
% %%Model 1
% model_type = 16;
% 
% par_free = [8.420147593529158 1.655254989751193e+02 round(1.202636718750000*180)];
% par_free = [8.420147593529158 1.655254989751193 round(1.202636718750000*180)];
% 
% call_optim_results(model_type, par_free, Ca_choice);
% 


% %%Model 2 %Model 242: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep, num_ves]
% model_type = 242;
% 
% par_free = [4.516909692627030 2 1.5e-6 4.098627212518727e+07 5.102599825769611e+02 1.592958707569855e+02 round(1.172753906249999*180)];
% 
% call_optim_results(model_type, par_free, Ca_choice);
% 


%%Model 3 %Model 42: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]

% model_type = 42;

% model_type = 46;
% 
% par_free = [13.7478 5 5.4999e-08 134.6764 236.5534 round(1.0038*180)];
% 
% call_optim_results(model_type, par_free, Ca_choice);
% 


% %Model S1 %Model 33: par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep, num_ves]
% model_type = 33;
% 
% par_free = [12.587204133650637 5 0.092900455549552 1.947727198548635e+02 10.697511762243380 1.412098817171219e+02 round(1.051855468750000*180)];
% 
% call_optim_results(model_type, par_free, Ca_choice);






% %Calcium simulation

% call_calcium_simulation;

