function [Ca_R_final, Ca_time, dists, Ca_80] = read_and_reshape_cylinder(filename, num_bins, rand_ves_on_off, par)

% Reads data from 3-column CalC file generated with
% plot 1D command. col1 = time, col2 = position, col3 = Ca conc.

%% Draw vesicles from new (gamma-ish) distribution

dists = determ_vesicle_distances(rand_ves_on_off, num_bins, par);

num_bins = length(dists);


%% Read and reshape CalC-file
CalC_plot = dlmread(filename,'');
time_1D = CalC_plot(:,1);
Position_1D = CalC_plot(:,2);
Ca_1D = CalC_plot(:,3);
Ca_1D = Ca_1D/1e6;
reshaped_Ca1D = reshape(Ca_1D, length(unique(Position_1D)), length(Ca_1D)/length(unique(Position_1D)));



num_grid = sum(time_1D(:,1)==time_1D(1,1));
times = time_1D(1:num_grid:length(time_1D),1);


[times_unique unique_inds] = unique(times);

reshaped_Ca1D_unique = reshaped_Ca1D(:,unique_inds);

Positions = unique(Position_1D);

Ca_R_vesicles = interp1(Positions(:,1),reshaped_Ca1D_unique,dists, 'PCHIP');

Ca_80_first = interp1(Positions(:,1),reshaped_Ca1D_unique,0.08,'PCHIP');
Ca_R_vesicles = Ca_R_vesicles';
Ca_times = times_unique/1000;
Ca_time = [0; Ca_times];
CaRest = Ca_R_vesicles(1,1);
Ca_R_final = ones(length(Ca_time), num_bins)*CaRest;
Ca_R_final(2:end,:) = Ca_R_vesicles;
Ca_80 = ones(length(Ca_time), 1)*CaRest;
Ca_80(2:end,:) = Ca_80_first;

end