function [cost_value] = determine_cost_peak1(peak1s)

load('DataSummary_extr.mat')

cost_value = sum(((peak1s - peak1_exp_mean).^2)./abs(peak1_exp_mean));
