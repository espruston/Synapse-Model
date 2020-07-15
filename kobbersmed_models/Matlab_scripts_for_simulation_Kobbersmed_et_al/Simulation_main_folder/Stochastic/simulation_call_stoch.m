function [peak1_means, peak2_extr_means, peak1_vars, peak2_extr_vars, ppr_extr_means, ttp_means, pVr1_means, stoch_EPSCs_cell, stoch_EPSC_means, peaksum_means, peaksum_vars, time_vector] =  simulation_call_stoch(par, CaExtracellular, rand_ves_on_off, save_data, savefilename, save_calc_loc, num_stim, stim_freq)

num_calc = length(CaExtracellular);
num_ves = par(22);
CalC_geometry = par(38);
AZ_size = par(16);
Q_max = par(17);
grid_points = par(34);
% num_stim = par(35);
height = par(36);
min_step = par(37);
num_sim = par(45);


if save_data == 3
    num_sim = 20
    par(45) = 20
elseif save_data == 66
    num_sim = 100;
end


stoch_EPSCs_cell = deal(cell(num_calc, num_sim));
[peak1s_all, peak2s_all, peaksums, pprs, ttps, pVr1s] = deal(zeros(num_calc, num_sim));

% act_states_cell = cell(num_calc, num_sim);
% ves_states_cell = cell(num_calc, num_sim);

Ca_residuals = zeros(num_calc, 1);
peak_extr_covs = zeros(num_calc, 1);

[Ca_R_equi_cell, Ca_time_equi_cell, time_vector_cell] = deal(cell(num_sim, 1));





if save_calc_loc == 1
    folder_name = './Sim_data/CalC_files/';
elseif save_calc_loc == 2
    folder_name = './Sim_data/new_calcium/';
end

for l = 1:num_calc

    filename = [folder_name 'Calcium_height' num2str(height) '_Qmax' num2str(Q_max) '_Ca' num2str(CaExtracellular(l)) '_numstim' num2str(num_stim) '_freq' num2str(stim_freq) '_AZ' num2str(AZ_size*1e3) '_geo' num2str(CalC_geometry) '_grid' num2str(grid_points)];
    

    if CalC_geometry == 1
        [~, Ca_residual, ~, ~] = read_and_reshape_residual_cylinder(filename);
        [Ca_R_vesicles_quant, Ca_time_vesicles_quant, dists_vesicles_equi, ~] = read_and_reshape_cylinder(filename, num_ves, 5, par);
    elseif CalC_geometry == 2
        [~, Ca_residual] = read_and_reshape_residual_BOX(filename, grid_points);
        [Ca_R_vesicles_quant, Ca_time_vesicles_quant, dists_vesicles_equi, ~] = read_and_reshape_BOX(filename, num_ves, 5, grid_points, par);        
    end

    sample_length = round(Ca_time_vesicles_quant(end)/min_step);
    if l == 1
        stoch_EPSCs = zeros(num_sim, sample_length+1);
        stoch_EPSC_means = zeros(num_calc, sample_length+1);
    end

    
    parfor z = 1:num_sim

        if CalC_geometry == 1
            [Ca_R_vesicles, Ca_time_vesicles, dists_vesicles, ~] = read_and_reshape_cylinder(filename, num_ves, rand_ves_on_off, par);
        elseif CalC_geometry == 2  
            [Ca_R_vesicles, Ca_time_vesicles, dists_vesicles, ~] = read_and_reshape_BOX(filename, num_ves, rand_ves_on_off, grid_points, par)
        else
            [Ca_R_vesicles, Ca_time_vesicles, dists_vesicles] = deal(NaN)
        end
        
        [peak1, peak2, time_vector, stoch_EPSC, ppr, peaksum, ttp1,  peak2_extr, ppr_extr, pVr1, ~, ~, ~, ~, ~, ~, peaksum_extr, ~] = main_exocytosis_stoch(Ca_time_vesicles, Ca_R_vesicles, par, stim_freq);
        
        
        peak1s_all(l,z) = peak1*1e9;
        peak2s_all(l,z) = peak2*1e9;
        ppr_extrs_all(l,z) = ppr_extr;
        peak2_extrs_all(l,z) = peak2_extr*1e9;
        peaksums(l,z) = peaksum*1e9;
        peaksums_extr(l,z) = peaksum_extr*1e9;
        pprs(l,z) = ppr;
        ttps(l,z) = ttp1;
        pVr1s(l,z) = pVr1;
        stoch_EPSCs_cell{l,z} = stoch_EPSC*1e9;
        
%         if save_data == 3
%             stoch_time_cell{l,z} = stoch_time;
%             ves_states_cell{l,z} = ves_states;
%             act_states_cell{l,z} = act_states_summed;
%             total_fast_states_cell{l,z} = total_fast_states;
%             total_slow_states_cell{l,z} = total_slow_states;
%             fast_states_cell{l,z} = fast_states;
%             slow_states_cell{l,z} = slow_states;
%             stoch_EPSC_chang = [stoch_EPSC zeros(1, 1e6-length(stoch_EPSC))];
%         end

        stoch_EPSCs(z,:) = stoch_EPSC*1e9; 
        time_vector_cell{z} = time_vector
    end
    
    
    stoch_EPSC_test = stoch_EPSCs_cell{1,1};
    sample_length = length(stoch_EPSC_test);
    stoch_EPSCs = stoch_EPSCs(:, 1:sample_length);
    stoch_EPSC_means(l, :) = mean(stoch_EPSCs,1);
    Ca_residuals(l) = Ca_residual*1e9;
    Ca_R_equi_cell{l} = Ca_R_vesicles_quant*1e6;
    Ca_time_equi_cell{l} = Ca_time_vesicles_quant*1e3;
    
    covmatrix_extr = cov(peak1s_all(l,:), peak2_extrs_all(l,:));
    peak_extr_covs(l) = covmatrix_extr(1,2);
end



time_vector = time_vector_cell{1};
stoch_EPSC_means = stoch_EPSC_means(:, 1:sample_length);

peak1_means = mean(peak1s_all, 2);
peak2_NOTEXTR_means = mean(peak2s_all, 2);
peak2_extr_means = mean(peak2_extrs_all, 2);
peaksum_means = mean(peaksums, 2);
peaksum_extr_means = mean(peaksums_extr, 2);
ppr_NOTEXTR_means = nanmean(pprs, 2);
ppr_extr_means = nanmean(ppr_extrs_all,2);
ttp_means = mean(ttps, 2);
pVr1_means = mean(pVr1s, 2);
peak1_vars = var(peak1s_all, 0, 2);
peak2_vars = var(peak2s_all, 0, 2);
peak2_extr_vars = var(peak2_extrs_all, 0, 2);
ppr_extr_vars = nanvar(ppr_extrs_all, 0, 2);
ppr_vars = nanvar(pprs, 0, 2);
peaksum_vars = var(peaksums, 0, 2);
peaksum_extr_vars = var(peaksums_extr, 0, 2);


save_folder = './Sim_data/new_results/';
[savename] = generate_savename(savefilename, 1, rand_ves_on_off, save_data, num_stim, stim_freq)
result_filename = [save_folder savename];

if ((save_data == 1) || (save_data == 3))
    save(result_filename, 'peak1s_all', 'peak2_extrs_all', 'ppr_extrs_all', 'peak1_means', 'peak2_extr_means', 'ppr_extr_means', 'ttp_means', 'pVr1_means', 'peak1_vars', 'peak2_extr_vars', 'stoch_EPSCs_cell', 'stoch_EPSC_means','time_vector', 'Ca_R_equi_cell', 'Ca_time_equi_cell', 'Ca_residuals', 'height', 'dists_vesicles_equi', 'num_sim', 'rand_ves_on_off', 'par', 'peak_extr_covs', 'ppr_extr_vars', 'peaksum_extr_means', 'peaksum_extr_vars', 'ppr_NOTEXTR_means', 'peak2_NOTEXTR_means', '-v7.3')
elseif save_data == 2
    save(result_filename, 'peak1_means', 'peak2_extr_means', 'peaksum_means', 'ppr_extr_vars', 'ppr_extr_means', 'ttp_means', 'pVr1_means', 'peak1_vars', 'peak2_extr_vars', 'peaksum_vars', 'peaksum_extr_means', 'peaksum_extr_vars', 'stoch_EPSC_means','time_vector', 'Ca_residuals', 'height', 'num_sim', 'rand_ves_on_off', 'par',  '-v7.3')
elseif save_data == 999
    save(result_filename, 'peak1s_all', 'peak1_means', 'peak2_extr_means', 'peaksum_means', 'ppr_extr_means', 'ttp_means', 'pVr1_means', 'peak1_vars', 'peak2_extr_vars', 'peaksum_vars', 'stoch_EPSCs_cell', 'stoch_EPSC_means','time_vector', 'Ca_R_equi_cell', 'Ca_time_equi_cell', 'Ca_residuals', 'height', 'dists_vesicles_equi', 'num_sim', 'rand_ves_on_off', 'par', 'act_states_cell', 'peak_extr_covs', 'ppr_extr_vars', 'ppr_vars', 'peaksum_extr_means', 'peaksum_extr_vars', 'fast_states_cell', 'slow_states_cell', 'total_fast_states_cell', 'total_slow_states_cell', 'stoch_time_cell', 'ves_states_cell', '-v7.3')
elseif save_data == 66 || save_data == 666
    save(result_filename, 'peak1s_all', 'peak2_extrs_all', 'ppr_extrs_all', 'rand_ves_on_off', 'num_sim', 'par', '-v7.3')
end
