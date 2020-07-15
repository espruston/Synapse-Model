function call_simulation_mod2_fig52_therest()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k_d_seconds = [1.3 1.4 1.6 1.7 1.8 1.9 2]*1e-6;

for k_d_second = k_d_seconds

% num_ves_factor = 1.5739;
SS_coop = 2;

rand_ves_on_off = 1;
% stoch_on_off = 2;

pVr2_hack = 0;

if k_d_second == 1.3e-9
    Q_maxes = [3 4 5 6 7 8];
else
    Q_maxes = [2 3 4 5 6 7 8]; %[2 3 4 5 6 7 8]; %Leaving out Q_max = 3 for now
end
% Q_maxes = fliplr(Q_maxes); %Start with highest Q_max for better time estimation

save_data = 666;

s_values = [1 50 100 200 300:50:550 600 650 700 750 800]; %[1 50 100 200 300:50:550 
% s_values = fliplr(s_values);
k_on_values = [5e7 2.5e7 1e7 7.5e7 5e6 2.5e6 1e6 5e5];






% k_d_values = [1 1.2 1.4 1.6 1.8 2]*1e-6;

k_rep = 157.5804; %k_rep from best fit
CaExtracellular = [0.75 1.5 3 6 10];
% CaExt = 0.75;
save_calc_loc = 1;
% s_values = 305:5:500;
% num_iterations = 2;

% num_ves_factors_peak1 = zeros(length(Q_maxes), length(k_on_values), length(s_values));
% num_ves_factors_bothpeaks = zeros(length(Q_maxes), length(k_on_values), length(s_values));

sim_num_total = length(Q_maxes)*length(k_on_values)*length(s_values);

% load('Fig52_num_ves_factors.mat')
for l = 1:length(Q_maxes)
    CalC_on_off = 1;
    Q_max = Q_maxes(l);
        
        
        
%     for n = 1:length(k_d_values)
%         k_d_second = k_d_values(n);




        for z = 1:length(k_on_values)
            k_on = k_on_values(z);
            for k = 1:length(s_values)
                s_val = s_values(k);

                tic;

                sim_num = (l-1)*length(s_values)*length(k_on_values) + (z-1)*length(s_values) + k;

                model_type = 241;

                par_free = [Q_max SS_coop k_d_second k_on s_val k_rep];
                pars = [Q_max k_d_second k_on s_val k_rep];
                [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);

                disp(['CalC simulation.'])

                if (z == 1) && (k == 1)
                    CalC_on_off = 1;
                    testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
                end
                disp('Calc simulation done. Time: ')

                CalC_on_off = 0;
                disp(' ')
                disp(['Now running simulation #' num2str(sim_num) ' of ' num2str(sim_num_total)])

%                 if k_d_second == 1.5e-6
                    [~, ~, num_ves_factor_peak1, num_ves_factor_bothpeaks, ~, ~, ~, ~, ~] = cost_model2_Qmax_kon_s_corrected_numves_krep_fig52(pars, 0);
%                 else
%                     num_ves_factor_peak1 = 1;
%                 end
                
                par_free(end+1) = round(180*num_ves_factor_peak1);
                model_type = 242;
                stoch_on_off = 1;
                [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);



    %             if sum(s_val == [600 650 700 750 800])>0
    %                 kkk = 1;
    %             else
    %                 kkk = 0;
    %             end
    %             


    %             for kk = kkk:length(CaExtracellular)
                    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular(1), save_data, savefilename, save_calc_loc, pVr2_hack);
                    disp('Time: ')
                    toc
                    CalC_on_off = 0;
    %             end
                            disp(' ')

    %             num_ves_factors_peak1(l,z,k) = num_ves_factor_peak1;
    %             num_ves_factors_bothpeaks(l,z,k) = num_ves_factor_bothpeaks;
            end
        end
%     end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % num_ves_factor = 1.5739;
% SS_coop = 2;
% 
% rand_ves_on_off = 1;
% % stoch_on_off = 2;
% k_d_second = 1.5e-6;
% 
% 
% Q_maxes = [2 3 4 5 6 7 8]; %[2 3 4 5 6 7 8]; %Leaving out Q_max = 3 for now
% % Q_maxes = fliplr(Q_maxes); %Start with highest Q_max for better time estimation
% 
% save_data = 666;
% 
% s_values = [50 100 200 300:50:550 600 650 700 750 800]; %[1 50 100 200 300:50:550 
% % s_values = fliplr(s_values);
% k_on_values = [5e7 2.5e7 1e7 7.5e7 5e6 2.5e6 1e6 5e5];
% 
% k_d_values = [1 1.2 1.4 1.6 1.8 2]*1e-6;
% 
% k_rep = 157.5804; %k_rep from best fit
% CaExtracellular = [0.75 1.5 3 6 10];
% % CaExt = 0.75;
% save_calc_loc = 1;
% % s_values = 305:5:500;
% % num_iterations = 2;
% 
% % num_ves_factors_peak1 = zeros(length(Q_maxes), length(k_on_values), length(s_values));
% % num_ves_factors_bothpeaks = zeros(length(Q_maxes), length(k_on_values), length(s_values));
% 
% sim_num_total = length(Q_maxes)*length(k_on_values)*length(s_values) * length(k_d_values);
% 
% % load('Fig52_num_ves_factors.mat')
% for l = 1:length(Q_maxes)
%     CalC_on_off = 1;
%     Q_max = Q_maxes(l);
%         
%         
%         
%     for n = 1:length(k_d_values)
%         k_d_second = k_d_values(n);
% 
% 
% 
% 
%         for z = 1:length(k_on_values)
%             k_on = k_on_values(z);
%             for k = 1:length(s_values)
%                 s_val = s_values(k);
% 
%                 tic;
% 
%                 sim_num = (l-1)*length(s_values)*length(k_on_values)*length(k_d_values) + (n-1)*length(s_values)*length(k_on_values) +  (z-1)*length(s_values) + k;
% 
%                 model_type = 241;
% 
%                 par_free = [Q_max SS_coop k_d_second k_on s_val k_rep];
%                 pars = [Q_max k_d_second k_on s_val k_rep];
%                 [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);
% 
%                 disp(['CalC simulation.'])
% 
%                 if (z == 1) && (k == 1) && (n == 1)
%                     CalC_on_off = 1;
%                     testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack);
%                 end
%                 disp('Calc simulation done. Time: ')
% 
%                 CalC_on_off = 0;
%                 disp(' ')
%                 disp(['Now running simulation #' num2str(sim_num) ' of ' num2str(sim_num_total)])
% 
%                 [~, ~, num_ves_factor_peak1, num_ves_factor_bothpeaks, ~, ~, ~, ~, ~] = cost_model2_Qmax_kon_s_corrected_numves_krep_fig52(pars, 0);
% 
%                 par_free(end+1) = round(180*num_ves_factor_peak1);
%                 model_type = 242;
%                 stoch_on_off = 1;
%                 [par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);
% 
% 
% 
%     %             if sum(s_val == [600 650 700 750 800])>0
%     %                 kkk = 1;
%     %             else
%     %                 kkk = 0;
%     %             end
%     %             
% 
% 
%     %             for kk = kkk:length(CaExtracellular)
%                     testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular(1), save_data, savefilename, save_calc_loc, pVr2_hack);
%                     disp('Time: ')
%                     toc
%                     CalC_on_off = 0;
%     %             end
%                             disp(' ')
% 
%     %             num_ves_factors_peak1(l,z,k) = num_ves_factor_peak1;
%     %             num_ves_factors_bothpeaks(l,z,k) = num_ves_factor_bothpeaks;
%             end
%         end
%     end
% end
% % save('Fig52_num_ves_factors.mat', 'num_ves_factors_peak1', 'num_ves_factors_bothpeaks')
% 
end