function [act_rate, inact_rate] = determine_activation_rates(act_rate_const, inact_rate_const, calcium, act_model_type)



    act_rate = calcium.^act_model_type * act_rate_const;
    inact_rate = inact_rate_const;
