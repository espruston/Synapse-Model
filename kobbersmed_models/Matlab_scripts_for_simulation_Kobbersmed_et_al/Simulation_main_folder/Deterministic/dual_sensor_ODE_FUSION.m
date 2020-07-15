function dual_SS = dual_sensor_ODE_FUSION(act_model_type, SS_PM, Ca_prim_type, Ca_R_sim, dist, input_pop, par, Ca_time_vesicles, t)

%Defines the rate equations solved with the ODE solver in deterministic
%simulations.


n_max = par(13);
m_max = par(14);

unprim_onestate = par(47);

if dist >= par(11)
    par(10) = 0;
end


num_states = (n_max + 1) * (m_max+1);
num_states_total = 3*(n_max+1)*(m_max+1);

if SS_PM == 0
    empty_index1 = num_states_total + 1;
    empty_index2 = num_states_total + 2;
    empty_index3 = num_states_total + 3 ;
    fused_index = num_states_total + 4;
elseif SS_PM == 1
    empty_index1 = (1:(m_max+1)) + num_states_total;
    empty_index2 = (((m_max+1)+1):2*(m_max+1)) + num_states_total ;
    empty_index3 = ((2*(m_max+1)+1):3*(m_max+1)) + num_states_total;
    fused_index = 3*(m_max+1)+1 + num_states_total;      
end

num_extra_states = length([empty_index1 empty_index2 empty_index3 fused_index]);

dual_SS = zeros(num_states_total+num_extra_states, 1);

Calcium = interp1(Ca_time_vesicles, Ca_R_sim, t);


input_pop_active_vector = input_pop(1:num_states);
input_pop_delay_vector = input_pop((num_states+1):2*num_states);
input_pop_inactive_vector = input_pop((2*num_states+1):3*num_states);

input_pop_active = reshape(input_pop_active_vector, m_max+1, n_max+1);
input_pop_delay = reshape(input_pop_delay_vector, m_max+1, n_max+1);
input_pop_inactive = reshape(input_pop_inactive_vector, m_max+1, n_max+1);

%Read activation and replenishment rates
k_rep = par(10);
act_rate_const = par(26);
inact_rate_const = par(27);
u_val = par(43);
prim_rate_model5 = (u_val~=1)*par(32);


[act_rate, inact_rate] = determine_activation_rates(act_rate_const, inact_rate_const, Calcium, act_model_type); %alpha and beta



delay_rate = par(28);
invdelay_rate = par(29);

%%Add calcium-independent replenishment
rep_rate_active = (k_rep + prim_rate_model5) * input_pop(empty_index1);
rep_rate_delay = (k_rep + prim_rate_model5) * input_pop(empty_index2);
rep_rate_inactive = (k_rep + prim_rate_model5) * input_pop(empty_index3);

[dual_SS_active_reshaped, fusion_rate] = calculate_propensities(Ca_R_sim, input_pop_active, par, Ca_time_vesicles, t, 1, SS_PM);
[dual_SS_delay_reshaped, ~] = calculate_propensities(Ca_R_sim, input_pop_delay, par, Ca_time_vesicles, t, 0, SS_PM);
[dual_SS_inactive_reshaped, ~] = calculate_propensities(Ca_R_sim, input_pop_inactive, par, Ca_time_vesicles, t, 0, SS_PM);


%%Add replenishment rate. Length of empty_index depends on SS_PM
dual_SS_active_reshaped(1:length(empty_index1),1) = dual_SS_active_reshaped(1:length(empty_index1),1) + rep_rate_active; %Adding replenishment rate
dual_SS_delay_reshaped(1:length(empty_index2),1) = dual_SS_delay_reshaped(1:length(empty_index2),1) + rep_rate_delay;
dual_SS_inactive_reshaped(1:length(empty_index3),1) = dual_SS_inactive_reshaped(1:length(empty_index3),1) + rep_rate_inactive;


dual_SS_active_noact = reshape(dual_SS_active_reshaped, 1, num_states)';
dual_SS_delay_noact = reshape(dual_SS_delay_reshaped, 1, num_states)';
dual_SS_inactive_noact = reshape(dual_SS_inactive_reshaped, 1, num_states)';
dual_SS_noact = [dual_SS_active_noact; dual_SS_delay_noact; dual_SS_inactive_noact];
dual_SS_noact(empty_index1) = fusion_rate-rep_rate_active; %Rep. and fusion at active states
dual_SS_noact(empty_index2) = -rep_rate_delay; %Replenishment at delay states
dual_SS_noact(empty_index3) = -rep_rate_inactive; %Replenishment at inactive states
dual_SS(fused_index) = sum(fusion_rate);

%Add activation and delay
dual_SS(1:num_states) = dual_SS_noact(1:num_states) ...
                        + input_pop((num_states+1):2*num_states)*delay_rate ...
                        - input_pop(1:num_states)*invdelay_rate;
                    
dual_SS((num_states+1):2*num_states) = dual_SS_noact((num_states+1):2*num_states) ...
                                       - input_pop((num_states+1):2*num_states)*delay_rate ...
                                       + input_pop(1:num_states)*invdelay_rate ...
                                       - input_pop((num_states+1):2*num_states)*inact_rate ...
                                       + input_pop((2*num_states+1):3*num_states)*act_rate;

dual_SS((2*num_states+1):3*num_states) = dual_SS((2*num_states+1):3*num_states) ...
                                         - input_pop((2*num_states+1):3*num_states) * act_rate ...
                                         + input_pop((num_states+1):2*num_states)*inact_rate;
                                   
dual_SS(empty_index1) = dual_SS_noact(empty_index1) + input_pop(empty_index2) * delay_rate...
                                              - input_pop(empty_index1) * invdelay_rate;
                                                        
dual_SS(empty_index2) = dual_SS_noact(empty_index2) - input_pop(empty_index2) * delay_rate ...
                                              + input_pop(empty_index1) * invdelay_rate ...
                                              - input_pop(empty_index2) * inact_rate ...
                                              + input_pop(empty_index3) * act_rate;
                                                            
dual_SS(empty_index3) = dual_SS_noact(empty_index3) - input_pop(empty_index3) * act_rate ...
                                              + input_pop(empty_index2) * inact_rate;
                                          
                                          
%%Add SS binding to empty states
if SS_PM == 1
    empty1_SS_prop = calculate_SS_propensities(Ca_R_sim, input_pop(empty_index1), par, Ca_time_vesicles, t);
    empty2_SS_prop = calculate_SS_propensities(Ca_R_sim, input_pop(empty_index2), par, Ca_time_vesicles, t);
    empty3_SS_prop = calculate_SS_propensities(Ca_R_sim, input_pop(empty_index3), par, Ca_time_vesicles, t);
    
    dual_SS(empty_index1) = dual_SS(empty_index1) + empty1_SS_prop;
    dual_SS(empty_index2) = dual_SS(empty_index2) + empty2_SS_prop;    
    dual_SS(empty_index3) = dual_SS(empty_index3) + empty3_SS_prop;
end

%%Add calcium-inhibited unpriming
if u_val ~= 1
    unprim_rate_const = par(33);
    
    for kk = 1:(m_max+1)
        inds = kk:(m_max+1):(((n_max+1)-1)*(m_max+1)+kk);
        unprim_rates = unprim_rate_const*u_val^(-(kk-1))*input_pop(inds);
        dual_SS(inds) = dual_SS(inds) - unprim_rates;
        dual_SS(empty_index1) = dual_SS(empty_index1) + sum(unprim_rates);
    end
    

end



    
%%Add Priming/unpriming rates
if Ca_prim_type > 0
    prim_kM = par(31);
    prim_rate_const = par(32);
    unprim_rate_const = par(33);
    unprim_rate_const_0 = par(44);
    
    empty_inds = [empty_index1 empty_index2 empty_index3];
    
    [Caprim_rate, Caunprim_rate] = calculate_Caprim_rate(Calcium, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);
    
    prim_prop = Caprim_rate * input_pop(empty_inds);
    unprim_prop_vesicles = Caunprim_rate * input_pop(1:num_states_total);
    
    if unprim_onestate
        unprim_prop_vesicles([2:num_states (num_states+2):(2*num_states) (2*num_states+2):(3*num_states)]) = 0;
    end
    
    if SS_PM == 0;
        unprim_prop_emp1 = sum(unprim_prop_vesicles(1:num_states));
        unprim_prop_emp2 = sum(unprim_prop_vesicles((num_states+1):(2*num_states)));
        unprim_prop_emp3 = sum(unprim_prop_vesicles((2*num_states+1):3*(num_states)));
    elseif SS_PM == 1; %%sum according to slow sensor status
        unprim_prop_vesicles_reshaped1 = reshape(unprim_prop_vesicles(1:num_states), m_max+1, n_max+1);
        unprim_prop_vesicles_reshaped2 = reshape(unprim_prop_vesicles((num_states+1):(2*num_states)), m_max+1, n_max+1);
        unprim_prop_vesicles_reshaped3 = reshape(unprim_prop_vesicles((2*num_states+1):3*(num_states)), m_max+1, n_max+1);
        unprim_prop_emp1 = sum(unprim_prop_vesicles_reshaped1, 2);
        unprim_prop_emp2 = sum(unprim_prop_vesicles_reshaped2, 2);
        unprim_prop_emp3 = sum(unprim_prop_vesicles_reshaped3, 2);
    end
    
    
    %unpriming
    dual_SS(1:num_states_total) = dual_SS(1:num_states_total) - unprim_prop_vesicles; %Adding unpriming rate to vesicles
    dual_SS(empty_inds) = dual_SS(empty_inds) + [unprim_prop_emp1; unprim_prop_emp2; unprim_prop_emp3]; %Adding unpriming rate to empty states
    %priming
    dual_SS(empty_inds) = dual_SS(empty_inds) - prim_prop; %

    Caprim_inds_first = 1:length(empty_index1);
    
%     if SS_PM == 0
%         Caprim_inds = 1:length([empty_index1 empty_index2 empty_index3]);
%     elseif SS_PM == 1
        Caprim_inds = [Caprim_inds_first (num_states + Caprim_inds_first) (2* num_states + Caprim_inds_first)];
%     end
        dual_SS(Caprim_inds) = dual_SS(Caprim_inds) + prim_prop;

end


