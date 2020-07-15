function run_more_reps(num_iterations, par_free, model_type, CalC_on_off, stoch_on_off, rand_ves_on_off, CaExtracellular, save_calc_loc, stim_freq, num_stim)

%Calls many repetitions of stochastic simulations in sets of 100. 
% Used for accurate result graphs. 
% Generates temporary result files by iteration, loads these and puts them together in one result file. 

%%Variables
%num_iterations: Number of sets of 100 simulations to be run.
%All other variables described 

save_data = 66;

num_sim = 100;

num_calc = length(CaExtracellular);
num_sim_total = num_sim*num_iterations;

[peak1s_all_coll, peak2_extrs_all_coll, ppr_extrs_all_coll] = deal(zeros(num_calc, num_sim_total));

save_folder = './Sim_data/new_results/';

pVr2_hack = 0;

for k = 1:num_iterations
    if k==1 && CalC_on_off == 1
        CalC_on_off_here = 1;
    else
        CalC_on_off_here = 0;
    end
    
    disp('Iteration: ')
    k
    
    
    [par_init, savefilename] = parameter_choices(par_free, model_type, k, rand_ves_on_off);


    tic;
    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off_here, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack, stim_freq, num_stim);
    toc
    
    savefilename = [savefilename '_Ca' num2str(CaExtracellular)];
    
    [savename] = generate_savename(savefilename, stoch_on_off, rand_ves_on_off, save_data);

    load([save_folder savename])
    
    peak1s_all_coll(:, ((k-1)*num_sim +1):k*num_sim) = peak1s_all;
    peak2_extrs_all_coll(:, ((k-1)*num_sim +1):k*num_sim) = peak2_extrs_all;
    ppr_extrs_all_coll(:, ((k-1)*num_sim +1):k*num_sim) = ppr_extrs_all;
    
    system(['rm ' save_folder savename]);
    
end

[~, savefilename_init] = parameter_choices(par_free, model_type, 0, rand_ves_on_off);



if length(CaExtracellular) == 1
    savefilename = [savefilename_init '_morereps' num2str(num_iterations) '_' num2str(num_sim) '_Ca' num2str(CaExtracellular) '_numstim' num2str(num_stim) '_stimfreq' num2str(stim_freq)]
else
    savefilename = [savefilename_init '_morereps' num2str(num_iterations) '_numstim' num2str(num_stim) '_stimfreq' num2str(stim_freq)]
end

save_name = generate_savename(savefilename, stoch_on_off, rand_ves_on_off, save_data);

if model_type == 46
    save_name = strrep(save_name, '.mat', '_MOD46.mat');
end


peak1_means = mean(peak1s_all_coll, 2);
peak2_extr_means = mean(peak2_extrs_all_coll, 2);
ppr_extr_means = nanmean(ppr_extrs_all_coll,2);
peak1_vars = var(peak1s_all_coll, 0, 2);
peak2_extr_vars = var(peak2_extrs_all_coll, 0, 2);
ppr_extr_vars = nanvar(ppr_extrs_all_coll, 0, 2);

save([save_folder save_name], 'ppr_extr_vars', 'peak2_extr_vars', 'peak1_vars', 'ppr_extr_means', 'peak1s_all_coll', 'peak2_extrs_all_coll', 'ppr_extrs_all_coll', 'peak1_means', 'peak2_extr_means', 'rand_ves_on_off', 'num_sim', 'par', '-v7.3')





