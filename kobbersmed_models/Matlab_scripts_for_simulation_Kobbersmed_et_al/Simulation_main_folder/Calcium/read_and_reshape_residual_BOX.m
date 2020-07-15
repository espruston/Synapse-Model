function [dist_residual, Ca_residual] = read_and_reshape_residual_BOX(filename, grid_points)

% Reads data from 3-column CalC file generated with
% plot 1D command. col1 = time, col2 = position, col3 = Ca conc.

%% Draw vesicles from rayleigh distribution, until all vesicles are within certain distance

num_grid = grid_points;

dist_residual = 100*1e-3;
dists_8ms = (1:370)*1e-3;

time_residual = 8*1e-3;




%% Read and reshape CalC-file
CalC_plot = dlmread(filename,'');
time_1D = CalC_plot(:,1);
Position_1D = CalC_plot(:,2);
Ca_1D = CalC_plot(:,3);
Ca_1D = Ca_1D/1e6;
reshaped_Ca1D = reshape(Ca_1D, length(unique(Position_1D)), length(Ca_1D)/length(unique(Position_1D)));


Times = time_1D(1:num_grid:length(time_1D),1);

Positions = unique(Position_1D);

Ca_R_interp = interp1(Positions(:,1),reshaped_Ca1D,dist_residual, 'PCHIP');

Ca_residual =  interp1(Ca_times,Ca_R_interp', time_residual, 'PCHIP');;
Times = Times/1000;

end