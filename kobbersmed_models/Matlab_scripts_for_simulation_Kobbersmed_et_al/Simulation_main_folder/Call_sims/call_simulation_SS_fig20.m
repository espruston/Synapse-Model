function call_simulation_SS_fig20()


% num_ves_factor = 1.5739;
SS_coop = 2;

rand_ves_on_off = 1;
stoch_on_off = 2;

% opt_Q_max = 5.3534;
Q_maxes = 3;

save_data = 2;

s_values = [1 5:5:500];
% s_values = 305:5:500;

for l = 1:length(Q_maxes)

    Q_max = Q_maxes(l)


    CalC_on_off = 1;
    disp(['CalC simulation.'])

    par_free = [Q_max SS_coop 1];
    model_type = 22;
    [par_init, savefilename] = parameter_choices(par_free, model_type);
    tic;
    testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
    disp('CalC simulation time:')
    toc


        CalC_on_off = 0;
    for k = 1:length(s_values)
        disp(['Simulation ' num2str(k) ' out of ' num2str(length(s_values)) '.'])
        opt_s = s_values(k);
        disp(['s value is ' num2str(opt_s) '.'])
        
        par_free = [Q_max SS_coop opt_s];
        model_type = 22;
        [par_init, savefilename] = parameter_choices(par_free, model_type);
        tic;
        testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
        disp('Exocytosis simulation time:')
        toc
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% % num_ves_factor = 1.5739;
% SS_coop = 6;
% 
% rand_ves_on_off = 1;
% stoch_on_off = 2;
% Q_max = 3;
% % opt_Q_max = 5.3534;
% 
% save_data = 2;
% 
% % s_values = [1 2.5 5:5:300];
% s_values = [2 3 4 5:5:40 50];
% CalC_on_off = 0;
% 
%     for k = 1:length(s_values)
%         disp(['Simulation ' num2str(k) ' out of ' num2str(length(s_values)) '.'])
%         opt_s = s_values(k);
%         disp(['s value is ' num2str(opt_s) '.'])
%         
%         par_free = [Q_max SS_coop opt_s];
%         model_type = 22;
%         [par_init, savefilename] = parameter_choices(par_free, model_type);
%         tic;
%         testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
%         disp('Exocytosis simulation time:')
%         toc
%     end


