function call_pVr_est_simulation()

tic;
CalC_on_off = 1;
pVr2_hack = 1;
save_data = 1;
stoch_on_off = 0;
rand_ves_on_off = 1000;

model_type = 241;
%  Q_maxes = [8.4201 13.7478];
CaExtracellular = [0.75 1.5 3 6 10];
save_calc_loc = 1;


% Model 242: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep, num_ves]. SS_PM = 1
Q_max = 4.5549;
SS_coop = 2;
k_D = 1.5e-6;
k_on = 39970832.4038;
s = 502.1079;
% k_rep = 157.5804;
% num_ves = 212;
    par_free(1) = Q_max;
    par_free(2) = SS_coop;
    par_free(3) = k_D;
    par_free(4) = k_on;
    par_free(5) = s;
    par_free(6) = 0;
%     par_free(7) = num_ves;


    [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);

    tic;
    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
    toc







%Model 1
par_free = [];

Q_max = 8.4201;
k_rep = 0;

model_type = 15;


par_free(1) = Q_max;
par_free(2) = k_rep

    [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);

    tic;
    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
    toc

    
%Model 41: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
par_free = [];
    CalC_on_off = 1
model_type = 41;
Q_max = 13.7478;
Ca_prim_type = 0%5;
prim_kM = 0%5.4999e-8;
prim_rate_constant = 0;%134.6764;
unprim_rate_const = 0;%236.5534;
% num_ves = 212;
    par_free(1) = Q_max;
    par_free(2) = Ca_prim_type;
    par_free(3) = prim_kM;
    par_free(4) = prim_rate_constant;
    par_free(5) = unprim_rate_const;
    

    %     par_free(7) = num_ves;
    [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);

    tic;
    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
    toc


%Model 32: par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep]
par_free = [];

model_type = 32;
Q_max = 12.5872;
act_model_type = 5;
beta = 0.0929;
gamma = 194.7727;
delta = 10.6975;
k_rep = 141.2099;
% num_ves = 212;
par_free(1) = Q_max;
par_free(2) = 0;
par_free(3) = 0;
par_free(4) = 0;
par_free(5) = 0;
par_free(6) = 0;

%     par_free(7) = num_ves;
[par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);

tic;
testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
toc
