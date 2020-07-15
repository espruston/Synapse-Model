function [peak1, peak2, time_vector, stoch_EPSC, ppr, peaksum, ttp2080, peak2_extr, ppr_extr, pVr1, stoch_time, act_states_summed, fast_states, slow_states, total_fast_states, total_slow_states, peaksum_extr, ves_states] = main_exocytosis_stoch(Ca_time_vesicles, Ca_R_vesicles, par, stim_freq)
%%%%%%%%%%%%%%%
%%%This script is the main exocytosis script. It contains the Gillespie
%%%algorithm call and the convolution of fused vesicles over time with the
%%%mEPSC. This is for paired pulse.
%%%%%%%%%%%%%%%

collect_states = 0

num_ves = par(22);
Ca_prim_type = par(39);
num_stim = par(35);
SS_PM = par(40);
act_model_type = par(42);
n_max = par(13);
m_max = par(14);
size_of_mini = par(19);
act_rate_const = par(26);
inact_rate_const = par(27);
delay_rate = par(28);
invdelay_rate = par(29);
u_val = par(43);
unprim_onestate = par(47);

Ca_bas = Ca_R_vesicles(1,1);
num_states_total = (n_max+1)*(m_max+1);

fuse_state = 100;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialise some simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sample_interval = 1e-6;
simulation_time = 0:sample_interval:Ca_time_vesicles(end);


tspan = [simulation_time(1) simulation_time(end)]; % Time window for the execution of directMethod


% Find steady state and vesicle states according to AZ_num random numbers:
if u_val == 1
    [SS_normal] = calculate_steady_state(par, Ca_bas, 1);
else
    SS_normal = zeros(num_states_total+1,1);
    steady_priming = calculate_steady_state_priming(par, Ca_bas); %steady_primin = [P, R0, R1, R2]; Total sum: 1.

    par_alt = par;
    par_alt([6:7 14]) = 0;
    SS_first_sensor = calculate_steady_state(par_alt, Ca_bas, 1); %First sensor steady
    for kk = 1:(m_max+1)
        inds = kk:(m_max+1):(((n_max+1)-1)*(m_max+1)+kk);
        SS_normal(inds) = SS_first_sensor'*steady_priming(kk+1);
    end
    SS_normal(num_states_total+1) = steady_priming(1);
end


if unprim_onestate
    prim_kM = par(31);
    prim_rate_const = par(32);
    unprim_rate_const = par(33);
    unprim_rate_const_0 = par(44);
    [Caprim_rate, Caunprim_rate] = calculate_Caprim_rate(Ca_bas, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);
    SS_normal(end+1) = SS_normal(1)*(Caunprim_rate/Caprim_rate);
    SS_normal = SS_normal/sum(SS_normal);
end

cumulative_steadyvesicles = cumsum(SS_normal);
ss_react_number = zeros(1, num_ves);

rnd_ss = rand(1,num_ves); 
for ff = 1:length(rnd_ss)
    ss_react_number(1,ff) = find(cumulative_steadyvesicles >= rnd_ss(1, ff),1, 'first');
end

% X is the vesicle state, n ranges from 0 to n_max for the fast sensor, and
% m ranges from 0 to m_max for the slow sensor.
X = ss_react_number;


n_ves = (X~=(num_states_total+1))'.*(ceil(X/(m_max+1)) -1)';
m_ves = (X~=(num_states_total+1))'.*(mod(X-1,m_max+1))';
ves_init_states = (10*m_ves + n_ves) + (X==(num_states_total+1))'*fuse_state;




%    ves_init_states = ves_init_states - ves_init_states; %NO STEADY STATE SOLUTION






%Add empty vesicles according to steady state steady state

if Ca_prim_type > 0 && ~unprim_onestate
    rnd_prim = rand(1,num_ves); 
    
    prim_kM = par(31);
    prim_rate_const = par(32);
    unprim_rate_const = par(33);
    unprim_rate_const_0 = par(44);
    [Caprim_rate, Caunprim_rate] = calculate_Caprim_rate(Ca_bas, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);
    prim_frac = Caprim_rate/(Caprim_rate + Caunprim_rate);
    unprim_frac = 1 - prim_frac;

    empt_ss = zeros(num_ves,1);
    for ee = 1:num_ves
        empt_ss(ee) = find([prim_frac 1] >= rnd_prim(ee), 1, 'first') - 1;
    end
    
    ves_init_states = ves_init_states - empt_ss.*n_ves + (SS_PM - 1) * m_ves + (empt_ss * 100);
    
end


%Determine initial activation



if (act_model_type ~= 0)
    [act_rate, inact_rate] = determine_activation_rates(act_rate_const, inact_rate_const, Ca_bas, act_model_type); %alpha and beta
    norm_const = 1 + (inact_rate / act_rate) + (delay_rate/invdelay_rate);

    inact_state = (inact_rate / act_rate) / norm_const;
    delay_state = 1 / norm_const;
    act_state = (delay_rate/invdelay_rate) / norm_const;

    cumulative_act_ss = cumsum([inact_state delay_state act_state]);

    ss_act_number = zeros(1, num_ves);
    rnd_act_ss = rand(1,num_ves); 
    for gg = 1:length(rnd_act_ss)
        ss_act_number(1,gg) = find(cumulative_act_ss >= rnd_act_ss(gg),1, 'first') - 2;
    end
else
    ss_act_number = ones(1, num_ves);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Stochastic exocytosis simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_vector = simulation_time';



[stoch_time, fused_ves, ves_states, act_states] = ...
    directMethod_newAlg_steps_SS(tspan, ves_init_states, ss_act_number, par, Ca_time_vesicles, Ca_R_vesicles, n_max, m_max, fuse_state, SS_PM, num_ves, act_model_type, Ca_prim_type, collect_states);


act_states_summed = zeros(3, length(stoch_time));
for v = 1:3
    act_states_summed(v,:) = sum(act_states == v-2, 1);
end

[fast_states, slow_states, total_fast_states, total_slow_states] = determine_total_slowfast_STOCH(ves_states, par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Convolving with mEPSC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fuse_times_appr = round(stoch_time(fused_ves > 0)*1e6); %In us, same index as time_vector

fuse_ves_vec = zeros(length(simulation_time), 1);

for j = 1:length(fuse_times_appr)
   fuse_ves_vec(fuse_times_appr(j)+1) = fuse_ves_vec(fuse_times_appr(j)+1) + 1;
end

pVr1 = sum(fuse_ves_vec(1:10000)) / num_ves;


[~, mEPSC_current] = smooth_mEPSC;

mEPSC_init = mEPSC_current;
mEPSC = -(mEPSC_init/min(mEPSC_init))*size_of_mini;


stoch_EPSC = conv(fuse_ves_vec, mEPSC); 
stoch_EPSC = stoch_EPSC(1:length(time_vector));

time_vector_size = size(time_vector);
EPSC_size = size(stoch_EPSC);

if time_vector_size ~= EPSC_size
    stoch_EPSC = stoch_EPSC';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if num_stim > 1
    [peak1, peak2, peaksum, ppr, peak2_extr, ppr_extr, ttp2080, peaksum_extr] = determine_peak_ppr(time_vector, stoch_EPSC, stim_freq);
else
    [peak2, peaksum, ppr, peak2_extr, ppr_extr, ttp2080] = deal(NaN);
end


warning on % Switch the warning statements on again