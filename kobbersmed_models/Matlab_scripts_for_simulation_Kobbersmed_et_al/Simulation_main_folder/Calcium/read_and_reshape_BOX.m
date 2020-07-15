function [Ca_R_vesicles, Times, dists, Ca_80] = read_and_reshape_BOX(filename, num_bins, rand_ves_on_off, grid_points, par)

% Reads data from 3-column CalC file generated with
% plot 1D command. col1 = time, col2 = position, col3 = Ca conc.

%% Draw vesicles from rayleigh distribution, until all vesicles are within certain distance

num_grid = grid_points;

dists = determ_vesicle_distances(rand_ves_on_off, num_bins, par);


%% Read and reshape CalC-file
CalC_plot = dlmread(filename,'');
time_1D = CalC_plot(:,1);
Position_1D = CalC_plot(:,2);
Ca_1D = CalC_plot(:,3);
Ca_1D = Ca_1D/1e6;
reshaped_Ca1D = reshape(Ca_1D, length(unique(Position_1D)), length(Ca_1D)/length(unique(Position_1D)));


Times = time_1D(1:num_grid:length(time_1D),1);

Positions = unique(Position_1D);

Ca_R_vesicles = interp1(Positions(:,1),reshaped_Ca1D,dists);

Ca_80 = interp1(Positions(:,1),reshaped_Ca1D,0.08);
Ca_R_vesicles = Ca_R_vesicles';
Times = Times/1000;

end