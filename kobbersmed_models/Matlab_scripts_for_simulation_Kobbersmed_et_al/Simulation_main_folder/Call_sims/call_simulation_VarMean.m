function call_simulation_VarMean()

%Model 11: par_free = Q_max
%Model 12: par_free = [Q_max kM]
%Model 13: par_free = [Q_max], k_rep = 0
%Model 14: par_free = [Q_max, num_ves], k_rep = 181.82
%Model 21: par_free = [Q_max, SS_coop, s]. SS_PM = 0.
%Model 22: par_free = [Q_max, SS_coop, s]. SS_PM = 1.
%Model 23: par_free = [Q_max, SS_coop, k_4, k_min4, s]. SS_PM = 0
%Model 24: par_free = [Q_max, SS_coop, k_4, k_min4, s]. SS_PM = 1
%Model 25: par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 0
%Model 26: par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 1
%Model 31: par_free = [Q_max, act_model_type, kMact, gamma, delta]
%Model 41: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
%Model 42: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
%Model 51: par_free = [Q_max, SS_coop, prim_rate, unprim_rate, u]; %Regehr values, SS_PM = 0.
%Model 51: par_free = [Q_max, SS_coop, k_4, k_min4, u]. SS_PM = 0


%%%%%%%

tic;
CalC_on_off = 1;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;


par_free(1) = 5.7001;
par_free(2) = round(180*1.6899);
model_type = 14;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc

%%%%%%%
clearvars;

tic;
CalC_on_off = 0;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;


par_free(1) = 5.7001;
par_free(2) = 2;
par_free(3) = 245;
par_free(4) = round(180*1.6899);
model_type = 26;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc

%%%%%%%
clearvars;

tic;
CalC_on_off = 0;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;


par_free(1) = 5.7001;
par_free(2) = 2;
par_free(3) = 245;
model_type = 22;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc

%%%%%%%
clearvars;

tic;
CalC_on_off = 1;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;

par_free(1) = 3;
par_free(2) = 2;
par_free(3) = 245;
model_type = 22;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc

%%%%%%%
clearvars;

tic;
CalC_on_off = 1;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;

par_free(1) = 4;
par_free(2) = 2;
par_free(3) = 245;
model_type = 22;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc

%%%%%%%
clearvars;

tic;
CalC_on_off = 1;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;

par_free(1) = 9.8578;
model_type = 11;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc

%%%%%%%
clearvars;

tic;
CalC_on_off = 1;
stoch_on_off = 1;
rand_ves_on_off = 1;

save_data = 2;

par_free(1) = 15;
par_free(2) = 5;
par_free(3) = 5e-8;
par_free(4) = 135;
par_free(5) = 375;
model_type = 41;


[par_init, savefilename] = parameter_choices(par_free, model_type);

testing_the_system_morereps(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);

disp('Simulation time: ')
toc


end

