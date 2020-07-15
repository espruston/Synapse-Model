function [dist_residual, Ca_residual, dists_10ms, Ca_10ms] = read_and_reshape_residual_cylinder(file_to_read)

% Reads data from 3-column CalC file generated with
% plot 1D command. col1 = time, col2 = position, col3 = Ca conc.

dist_residual = 122*1e-3;
dists_10ms = (1:370)*1e-3;

time_residual = 10*1e-3;


%% Read and reshape CalC-file
CalC_plot = dlmread(file_to_read,'');
time_1D = CalC_plot(:,1);
Position_1D = CalC_plot(:,2);
Ca_1D = CalC_plot(:,3);
Ca_1D = Ca_1D/1e6;
reshaped_Ca1D = reshape(Ca_1D, length(unique(Position_1D)), length(Ca_1D)/length(unique(Position_1D)));


num_grid = sum(time_1D(:,1)==time_1D(1,1));
times = time_1D(1:num_grid:length(time_1D),1);


[times_unique, unique_inds] = unique(times);

reshaped_Ca1D_unique = reshaped_Ca1D(:,unique_inds);

Positions = unique(Position_1D);

Ca_times = times_unique/1000;
Ca_R_interp = interp1(Positions,reshaped_Ca1D_unique,dist_residual, 'PCHIP');
Ca_R_interp_morepositions = interp1(Positions,reshaped_Ca1D_unique,dists_10ms, 'PCHIP');
Ca_residual = interp1(Ca_times,Ca_R_interp', time_residual, 'PCHIP');
Ca_10ms = interp1(Ca_times, Ca_R_interp_morepositions', time_residual, 'PCHIP');


end