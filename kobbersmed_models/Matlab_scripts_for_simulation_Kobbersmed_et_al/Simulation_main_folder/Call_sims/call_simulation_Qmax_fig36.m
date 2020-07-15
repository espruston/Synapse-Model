function call_simulation_Qmax_fig36()


rand_ves_on_off = 1;
save_data = 66;
num_iterations = 10;

Q_maxes = [0.1:0.05:0.45 0.5:0.4:14.9];
Q_maxes = fliplr(Q_maxes);  %Start with highest Q_max for better time estimation

% Q_maxes = 0.01
CaExt = 0.75;
save_calc_loc = 1;

for k = 1:length(Q_maxes)
        disp(['Simulation ' num2str(k) ' out of ' num2str(length(Q_maxes)) '.'])
        tic;
        CalC_on_off = 1;
        par_free(1) = Q_maxes(k);
        par_free(2) = 165.5255;
        par_free(3) = 200;
        model_type = 16;
        [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
        
        stoch_on_off = 0;
        testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExt, save_data, savefilename, save_calc_loc);
        stoch_on_off = 1;
        CalC_on_off = 0;
        run_more_reps(num_iterations, par_free, model_type, CalC_on_off, stoch_on_off, rand_ves_on_off, CaExt, save_calc_loc)

        toc
        
        CalC_on_off = 1;
        par_free(1) = Q_maxes(k);
        par_free(2) = 0;
        par_free(3) = 200;
        model_type = 16;
        [par_init, savefilename] = parameter_choices(par_free, model_type, 0);
        stoch_on_off = 0;
        testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExt, save_data, savefilename, save_calc_loc);
        stoch_on_off = 1;
        run_more_reps(num_iterations, par_free, model_type, CalC_on_off, stoch_on_off, rand_ves_on_off, CaExt, save_calc_loc)

        toc        

end


