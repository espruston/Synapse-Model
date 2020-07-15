function call_simulation_mod4_fig54()


% num_ves_factor = 1.5739;
prim_type = 5;
prim_rate = 135;

rand_ves_on_off = 1;
% stoch_on_off = 2;
k_M_prims = 3e-8:0.5e-8:8e-8;

% opt_Q_max = 5.3534;
Q_maxes = [6 7 8 9 10 11 12 13 14 15];

save_data = 66;

unprim_rate_values = 100:50:1000;

CaExtracellular = [0.75 1.5 3 6 10];
% CaExt = 0.75;
save_calc_loc = 1;
% s_values = 305:5:500;
% num_iterations = 2;

% num_ves_factors_peak1 = zeros(length(Q_maxes), length(k_on_values), length(s_values));
% num_ves_factors_bothpeaks = zeros(length(Q_maxes), length(k_on_values), length(s_values));

sim_num_total = length(Q_maxes)*length(k_M_prims)*length(unprim_rate_values);

% load('Fig52_num_ves_factors.mat')

for l = 1:length(Q_maxes)
    CalC_on_off = 1;
    Q_max = Q_maxes(l);


    for z = 1:length(k_M_prims)
        k_M_prim = k_M_prims(z);
        for k = 1:length(unprim_rate_values)
            unprim_rate = unprim_rate_values(k);

            tic;

            sim_num = (l-1)*length(unprim_rate_values)*length(k_M_prims) + (z-1)*length(unprim_rate_values) + k;

            model_type = 41;

            par_free = [Q_max prim_type k_M_prim prim_rate unprim_rate];
            pars = [Q_max k_M_prim prim_rate unprim_rate];
            [par_init, savefilename] = parameter_choices(par_free, model_type, 0);

            disp(['CalC simulation.'])

%             if ((l == 1) && (z == 1) && (k == 1))
%                 CalC_on_off = 1;
                testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc);
%                 CalC_on_off = 0;
%             end
            disp('Calc simulation done. Time: ')
%             toc
            CalC_on_off = 0;
            
            disp(' ')

%             tic;
            disp(['Now running simulation #' num2str(sim_num) ' of ' num2str(sim_num_total)])

            [~, ~, num_ves_factor_peak1, num_ves_factor_bothpeaks, ~, ~, ~, ~, ~] = cost_model4_fig54(pars, 0);

            par_free(end+1) = round(180*num_ves_factor_peak1);
            model_type = 42;
            stoch_on_off = 1;
            [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
            
            
%             if sum(s_val == [600 650 700 750 800])>0
%                 kkk = 1;
%             else
%                 kkk = 0;
%             end
%             
            
            
%             for kk = kkk:length(CaExtracellular)
                testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular(1), save_data, savefilename, save_calc_loc);
                disp('Time: ')
                toc
                CalC_on_off = 0;
%             end
                        disp(' ')

%             num_ves_factors_peak1(l,z,k) = num_ves_factor_peak1;
%             num_ves_factors_bothpeaks(l,z,k) = num_ves_factor_bothpeaks;
        end
    end
end

% save('Fig52_num_ves_factors.mat', 'num_ves_factors_peak1', 'num_ves_factors_bothpeaks')

