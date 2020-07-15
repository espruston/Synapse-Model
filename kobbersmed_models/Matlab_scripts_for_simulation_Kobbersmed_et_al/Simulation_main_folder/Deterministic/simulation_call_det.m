function [pVr_bins_cell_10ms, time_vector, EPSCs, peak1s, peak2s_extr, pprs_extr, states_all_cell, pVr1s, ttps, act_1s, act_2s, act_diffs] = simulation_call_det(par, CaExtracellular, num_bins, rand_ves_on_off, save_data, savefilename, save_calc_loc, pVr2_hack, num_stim, stim_freq)

num_calc = length(CaExtracellular);
AZ_size = par(16);
Q_max = par(17);
grid_points = par(34);
% num_stim = par(35);
height = par(36);

n_max = par(13);
m_max = par(14);

CalC_geometry = par(38);


%Initializing variables
[act_1s, act_2s, act_diffs, pprs_NOTEXTR, peak1s, peak2s_NOTEXTR, peak2s_extr, pprs_extr, peaksums_extr, pVr1s, ttps, Ca_residuals] = deal(zeros(1, length(CaExtracellular)));
[slow_states_cell, fast_states_cell, total_fast_states_cell, total_slow_states_cell, states_all_cell, Ca_R_vesicles_cell, Ca_time_vesicles_cell, ves_avail_cell] =  deal(cell(length(CaExtracellular),1));
pVr_bins_cell_10ms = zeros(length(CaExtracellular),num_bins);


%Where to read calcium files from
if save_calc_loc == 1
    folder_name = './Sim_data/CalC_files/';
elseif save_calc_loc == 2
    folder_name = './Sim_data/new_calcium/';
end

for l = 1:num_calc

    filename = [folder_name 'Calcium_height' num2str(height) '_Qmax' num2str(Q_max) '_Ca' num2str(CaExtracellular(l)) '_numstim' num2str(num_stim) '_freq' num2str(stim_freq) '_AZ' num2str(AZ_size*1e3) '_geo' num2str(CalC_geometry) '_grid' num2str(grid_points)];

    
    if rand_ves_on_off ==1 %If stochastic SV placement then use distribtion bins
        rand_ves_on_off_here = 0;
    else
        rand_ves_on_off_here = rand_ves_on_off;
    end
    
    if CalC_geometry == 1 %cylinder geometry
        [~, Ca_residual, ~, ~] = read_and_reshape_residual_cylinder(filename);
        [Ca_R_vesicles, Ca_time_vesicles, dists_vesicles, ~] = read_and_reshape_cylinder(filename, num_bins, rand_ves_on_off_here, par);
    elseif CalC_geometry == 2 %box geometry
        [~, Ca_residual] = read_and_reshape_residual_BOX(filename, grid_points);
        [Ca_R_vesicles, Ca_time_vesicles, dists_vesicles, ~] = read_and_reshape_BOX(filename, num_bins, rand_ves_on_off_here, grid_points, par);        
    end
    
    
    
    [pVr1_bins_10ms, time_vector, EPSC, peak1, peak2, ppr, states_all, pVr1, ttp1,  peak2_extr, ppr_extr, fast_states, slow_states, total_fast_states, total_slow_states, peaksum_extr, pVr1_bins_3_5ms, pVr1_bins_4ms, pVr1_bins_4_5ms, pVr1_bins_5ms, pVr2_bins_20ms, pVr2_bins_13_5ms, pVr2_bins_14ms, pVr2_bins_14_5ms, pVr2_bins_15ms, ves_avail] = determining_states_with_fusion(par, Ca_time_vesicles, Ca_R_vesicles, dists_vesicles, pVr2_hack, stim_freq); 

    if pVr2_hack == 0
        [activation, ~, ~] = determine_activation_from_simdata(time_vector, states_all, 0, n_max, m_max);
        act_1 = activation(1,:);
        act_2 = activation(1e4,:);
        act_diff = act_2 - act_1;
    
        ves_avail_cell{l} = ves_avail;
    else
        [activation, act_1, act_2, act_diff] = deal(NaN);
    end

    
    
    
    if l == 1
        EPSCs = zeros(num_calc, length(EPSC));
    end
    


    
     
    pprs_NOTEXTR(l) = ppr;
    peak1s(l) = peak1*1e9;
    peak2s_NOTEXTR(l) = peak2*1e9;
    pVr1s(l) = pVr1;
    ttps(l) = ttp1;
    EPSCs(l,:) = EPSC*1e9;
    fast_states_cell{l} = fast_states;
    slow_states_cell{l} = slow_states;
    total_fast_states_cell{l} = total_fast_states;
    total_slow_states_cell{l} = total_slow_states;
    states_all_cell{l} = states_all;
    Ca_residuals(l) = Ca_residual*1e9;
    pVr_bins_cell_10ms(l,:) = pVr1_bins_10ms;
    pVr_bins_cell_3_5ms(l,:) = pVr1_bins_3_5ms;
    pVr_bins_cell_4ms(l,:) = pVr1_bins_4ms;
    pVr_bins_cell_4_5ms(l,:) = pVr1_bins_4_5ms;
    pVr_bins_cell_5ms(l,:) = pVr1_bins_5ms;

    pVr2_bins_cell_20ms(l,:) = pVr2_bins_20ms;
    pVr2_bins_cell_13_5ms(l,:) = pVr2_bins_13_5ms;
    pVr2_bins_cell_14ms(l,:) = pVr2_bins_14ms;
    pVr2_bins_cell_14_5ms(l,:) = pVr2_bins_14_5ms;
    pVr2_bins_cell_15ms(l,:) = pVr2_bins_15ms;

    
    
    
    Ca_R_vesicles_cell{l} = Ca_R_vesicles*1e6;
    Ca_time_vesicles_cell{l} = Ca_time_vesicles*1e3;
    act_1s(l) = act_1;
    act_2s(l) = act_2;
    act_diffs(l) = act_diff;
    peak2s_extr(l) = peak2_extr*1e9;
    pprs_extr(l) = ppr_extr;
    peaksums_extr(l) = peaksum_extr*1e9;
    

end



save_folder = './Sim_data/new_results/';
[savename] = generate_savename(savefilename, 0, rand_ves_on_off, save_data, num_stim, stim_freq);
result_filename = [save_folder savename];

if ((save_data == 1) || (save_data == 3)) %Many results including EPSCs.
    save(result_filename, 'act_1s', 'act_2s', 'pVr_bins_cell_10ms', 'pVr_bins_cell_3_5ms', 'pVr_bins_cell_4ms', 'pVr_bins_cell_4_5ms', 'pVr_bins_cell_5ms', 'pVr2_bins_cell_20ms', 'pVr2_bins_cell_13_5ms', 'pVr2_bins_cell_14ms', 'pVr2_bins_cell_14_5ms', 'pVr2_bins_cell_15ms', 'time_vector', 'EPSCs', 'peak1s', 'fast_states_cell', 'pVr1s', 'ttps', 'states_all_cell', 'slow_states_cell', 'Ca_residuals', 'Ca_R_vesicles_cell', 'Ca_time_vesicles_cell', 'total_fast_states_cell', 'total_slow_states_cell', 'peak2s_extr', 'pprs_extr', 'peaksums_extr', 'dists_vesicles', 'ves_avail_cell', 'pprs_NOTEXTR', 'peak2s_NOTEXTR', 'par', '-v7.3')
elseif save_data == 2 %Fewer results
    save(result_filename, 'act_1s', 'act_2s', 'time_vector', 'peak1s', 'pVr1s', 'ttps', 'Ca_residuals', 'peak2s_extr', 'pprs_extr', 'peaksums_extr', 'par', '-v7.3')
elseif save_data == 3 %Only the essential results
    save(result_filename, 'act_1s', 'act_2s', 'time_vector', 'EPSCs', 'peak1s', 'Ca_R_vesicles_cell', 'Ca_time_vesicles_cell', 'peak2s_extr', 'pprs_extr', 'dists_vesicles',  'par', '-v7.3')
elseif save_data == 666 %Only the essential results
    save(result_filename, 'peak1s', 'peak2s_extr', 'pprs_extr', 'par', '-v7.3')
end
end
