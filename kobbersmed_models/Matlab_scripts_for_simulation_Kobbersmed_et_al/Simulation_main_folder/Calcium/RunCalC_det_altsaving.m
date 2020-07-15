function []=RunCalC_det_altsaving(CaExtracellular, par, CaRest, long_sim_on_off, num_stim, stim_freq) %Bm_conc_model1, ATP_conc_model1, fast_buffer_on_off, fast_buffer_params,

%This script modifies the parameter file for calcium simulation according to the chosen parameters in testing_the_system.m
%Calcium files and parameter are save in './Sim_data/new_calcium/'
%Simulation time is 163.5 ms if long_sim_on_off== 1. Else 20.5 ms.
    
CalC_geometry = par(38);
CaRest_uM = CaRest*1e6;
grid_points = par(34);
% num_stim = par(35);
height = par(36);
hill = par(41);


Q_max = par(17);
K_m = par(18);
AZ_size = par(16);
z_dist = par(30);


iChannel = calculate_ichan(Q_max, K_m, CaExtracellular, hill);


if long_sim_on_off
    PAR = fileread(['./Calcium/InitFile_SSREP_numstim' num2str(num_stim) '_longsim.par']);
else
    PAR = fileread(['./Calcium/InitFile_SSREP_numstim' num2str(num_stim) '.par']);
end


iChannel2 = sprintf('iChannel = %.6f %%', iChannel);
newPAR = strrep(PAR, 'iChannel = ', iChannel2);


if CalC_geometry == 1
    geometry_text = 'geometry = cylindrical';
    newPAR = strrep(newPAR, 'geometry', geometry_text);
    gridtext = ['grid ' num2str(grid_points) ' N'];
    newPAR = strrep(newPAR, 'grid', gridtext);
elseif CalC_geometry == 2
    volumetext = ['volume -' num2str(AZ_size) ' ' num2str(AZ_size) ' -' num2str(AZ_size) ' ' num2str(AZ_size) ' 0 height'];
    newPAR = strrep(newPAR, 'volume 0 AZ_size 0 height', volumetext);
    gridtext = ['grid ' num2str(grid_points) ' ' num2str(grid_points) ' ' num2str(grid_points)];
    newPAR = strrep(newPAR, 'grid', gridtext);
    boundtext = 'Ca.bc Noflux Noflux Noflux Noflux ...';
    newPAR = strrep(newPAR, 'Ca.bc Noflux Noflux ...', boundtext);
    sourcetext = 'Ca.source 0 0 0';
    newPAR = strrep(newPAR, 'Ca.source 0 0 ', sourcetext);
    newPAR = strrep(newPAR, 'geometry', ' ');

end



nodes = round(height/0.01);
nodes = nodes + (mod(nodes,2)==0);

nodes = sprintf('N = %.f %%', nodes);
newPAR = strrep(newPAR, 'N = ', nodes);


if CalC_geometry == 1
    iterationA1 = ['plot 1D.mute Ca r '  num2str(z_dist) ' ./Sim_data/new_calcium/Calcium_height' num2str(height) '_Qmax' num2str(Q_max) '_Ca' num2str(CaExtracellular) '_numstim' num2str(num_stim) '_freq' num2str(stim_freq) '_AZ' num2str(AZ_size*1e3) '_geo' num2str(CalC_geometry) '_grid' num2str(grid_points)];
    newPAR = strrep(newPAR, 'plot 1D.mute Ca r 0.01 CalciumA_1D', iterationA1);
elseif CalC_geometry == 2
    iterationA1 = ['plot 1D.mute Ca x 0 '  num2str(z_dist) ' ./Sim_data/new_calcium/Calcium_height' num2str(height) '_Qmax' num2str(Q_max) '_Ca' num2str(CaExtracellular) '_numstim' num2str(num_stim) '_freq' num2str(stim_freq) '_AZ' num2str(AZ_size*1e3) '_geo' num2str(CalC_geometry) '_grid' num2str(grid_points)];
    newPAR = strrep(newPAR, 'plot 1D.mute Ca r 0.01 CalciumA_1D', iterationA1);    
end
    
   

CaRest_str = sprintf('Ca.bgr = %.5f %', CaRest_uM);
newPAR = strrep(newPAR, 'Ca.bgr = 0.05', CaRest_str);

height_str = sprintf('height = %.5f %', height);
newPAR = strrep(newPAR, 'height = ', height_str);

AZ_size_str = sprintf('AZ_size = %.5f %', AZ_size);
newPAR = strrep(newPAR, 'AZ_size = ', AZ_size_str);

numstim_str = sprintf('numstim = %.5f %', num_stim);
newPAR = strrep(newPAR, 'numstim = ', numstim_str);

filestr1 = ['./Sim_data/new_calcium/newpar_height' num2str(height) '_Qmax' num2str(Q_max) '_Ca' num2str(CaExtracellular) '_numstim' num2str(num_stim) '_AZ' num2str(AZ_size*1e3) '_geo' num2str(CalC_geometry) '_grid' num2str(grid_points)];


newParameterfile = fopen(filestr1,'w');
fprintf(newParameterfile,'%s\n',newPAR);
fclose(newParameterfile);


fclose('all');

%run executable with content of name as input
if strcmp(computer, 'PCWIN64') %machine_ID == 'PCWIN64' %ispc==1
   [~,~]=system(['.\calc.exe ' filestr1]);
elseif strcmp(computer, 'MACI64') %machine_ID == 'MACI64' %ismac&islinux==1
 [~,~]=system(['./cmac686 ' filestr1]);
elseif strcmp(computer, 'GLNXA64') %machine_ID == 'GLNXA64' %isunix==1
     [~,~] = system(['./calc ' filestr1]);
else
   disp('Platform not supported')
end

end