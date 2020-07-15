function [act_prop] = propensities_newAlg_act(act_number_vesicles, par, Ca_R_temp, act_model_type)

%Reaction rates of activation

act_rate_const = par(26);
inact_rate_const = par(27);
delay_rate = par(28);
invdelay_rate = par(29);
[act_rate, inact_rate] = determine_activation_rates(act_rate_const, inact_rate_const, Ca_R_temp, act_model_type); %alpha and beta

act_prop = (act_number_vesicles == -1) .* act_rate + (act_number_vesicles == 0) .* (inact_rate + delay_rate) + (act_number_vesicles == 1) .* invdelay_rate;