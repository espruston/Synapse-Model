function [cost_value_peak1, num_ves_factor_peak1] = determine_cost_peak1_numves_factor(peak1s)

load('DataSummary_extr.mat')

cost_function = @(x)sum(((peak1s*x - peak1_mean_all_cells_exp).^2)./abs(peak1_mean_all_cells_exp));

[num_ves_factor_peak1, cost_value_peak1] = fminsearch(cost_function, 2);