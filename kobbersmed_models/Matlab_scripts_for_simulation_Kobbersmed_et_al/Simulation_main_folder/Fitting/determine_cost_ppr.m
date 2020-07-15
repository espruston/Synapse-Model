function [cost_value] = determine_cost_ppr(pprs)

load('NewDataDrosophila.mat')

cost_value = sum(((pprs - ppr_exp_mean).^2)./abs(ppr_exp_mean));
