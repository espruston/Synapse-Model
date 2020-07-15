function [cost_value_bothpeaks, num_ves_factor_bothpeaks] = determine_cost_bothpeaks_numves_factor(peak1s, peak2s)

load('DataSummary_extr.mat')

cost_function = @(x)(sum(((peak1s*x - peak1_mean_all_cells_exp).^2)./abs(peak1_mean_all_cells_exp))...
                    +sum(((peak2s*x - peak2_extr_mean_all_cells_exp).^2)./abs(peak2_extr_mean_all_cells_exp)));

[num_ves_factor_bothpeaks, cost_value_bothpeaks] = fminsearch(cost_function, 2);