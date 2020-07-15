function [savename] = generate_savename(savefilename, stoch_on_off, rand_ves_on_off, save_data, num_stim, stim_freq)

savename = [savefilename '_num_stim' num2str(num_stim) 'stim_freq' num2str(stim_freq)];

if stoch_on_off == 0
    savename = [savename '_DET'];
elseif stoch_on_off == 1
    savename = [savename '_STOCH'];
end

if rand_ves_on_off == 1000
    savename = [savename '_pVr'];
end

if save_data == 1
    savename = [savename '.mat'];
elseif save_data == 2
    savename = [savename '_fewresults.mat'];
elseif save_data == 3
    savename = [savename '_fewrepetitions.mat'];
elseif sum(save_data == [66 666])
    savename = [savename '_fewerresults.mat'];
elseif save_data == 999
    savename = [savename '_manyresults.mat'];
end

