function [pVr1_bins_10ms, time_vector, EPSC, peak1, peak2, ppr, states_all, pVr1, ttp2080,  peak2_extr, ppr_extr, fast_states, slow_states, total_fast_states, total_slow_states, peaksum_extr, pVr1_bins_3_5ms, pVr1_bins_4ms, pVr1_bins_4_5ms, pVr1_bins_5ms, pVr2_bins_20ms, pVr2_bins_13_5ms, pVr2_bins_14ms, pVr2_bins_14_5ms, pVr2_bins_15ms, ves_avail] = determining_states_with_fusion(par, Ca_time_vesicles, Ca_R_vesicles, dists_vesicles, pVr2_hack, stim_freq)

Ca_bas = Ca_R_vesicles(1,1);


num_ves = par(22);
SS_PM = par(40);
SS_coop = par(14);
num_stim = par(35);
act_model_type = par(42);
Ca_prim_type = par(39);
u_val = par(43);
unprim_onestate = par(47);

n_max = par(13);
m_max = par(14);
size_of_mini = par(19);

num_states_total = 3*(n_max+1)*(m_max+1);


num_bins = length(Ca_R_vesicles(1,:));
num_states = (n_max+1) * (m_max+1);

act_rate_const = par(26);
inact_rate_const = par(27);
delay_rate = par(28);
invdelay_rate = par(29);





simulation_time = Ca_time_vesicles(end); %s
sample_rate = 1e-6;


num_empty_fus = 3*(SS_PM*SS_coop+1) + 1; %3 empty states (A,D,I) times number of SS states, if SS_Pm==1. +1 for fusion state. 
fuse_state = num_states_total + num_empty_fus;

%%%%%



num_ves_in_steady = num_ves/num_bins;



%Initializing
time_vector = 0:sample_rate:simulation_time;
slow_states = zeros(length(time_vector), (m_max+1));
fast_states = zeros(length(time_vector), (n_max+1));
fused_states = zeros(length(time_vector), num_bins);
empty_states_1 = zeros(length(time_vector), num_bins);
empty_states_2 = zeros(length(time_vector), num_bins);
empty_states_3 = zeros(length(time_vector), num_bins);
allstatesinallbins_cell = cell(1, num_bins);
states_all = zeros(length(time_vector), 3*num_states+num_empty_fus);
fused_states_EPSC = zeros(length(time_vector),1);



    %Steady state from direct calculation
[SS_delaystate_nonorm] = calculate_steady_state(par, Ca_bas, num_ves_in_steady);


%Site activation: Not all AZs are activated at the beginning
if act_model_type == 0 %Not site activation model
    SS_activestate_nonorm = SS_delaystate_nonorm;
    SS_inactstate_nonorm = zeros(size(SS_delaystate_nonorm));
    SS_delaystate_nonorm = zeros(size(SS_inactstate_nonorm));
    par(26:29) = 0;
else %Site activation model
    [act_rate, inact_rate] = determine_activation_rates(act_rate_const, inact_rate_const, Ca_bas, act_model_type); %alpha and beta
    SS_inactstate_nonorm = (inact_rate / act_rate) * SS_delaystate_nonorm;
    SS_activestate_nonorm = (delay_rate/invdelay_rate) * SS_delaystate_nonorm;
end

    steadyvesicles_reshaped_nonorm = [SS_activestate_nonorm; SS_delaystate_nonorm; SS_inactstate_nonorm];
    steadyvesicles_reshaped = steadyvesicles_reshaped_nonorm * (num_ves_in_steady/sum(sum(steadyvesicles_reshaped_nonorm)));
    steadyvesicles_reshaped(end+(1:num_empty_fus)) = 0; %Empty site state 1:3 and fusion state
    
%Unpriming model
if Ca_prim_type > 0 || u_val ~= 1
    prim_kM = par(31);
    prim_rate_const = par(32);
    unprim_rate_const = par(33);
    unprim_rate_const_0 = par(44);

    if Ca_prim_type > 0
        [Caprim_rate, Caunprim_rate] = calculate_Caprim_rate(Ca_bas, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);
        prim_frac = Caprim_rate/(Caprim_rate + Caunprim_rate);
        unprim_frac = 1 - prim_frac;

        if SS_PM == 0
            SS_empty = [sum(steadyvesicles_reshaped(1:num_states)) sum(steadyvesicles_reshaped((num_states+1):2*num_states)) sum(steadyvesicles_reshaped((2*num_states+1):3*num_states))] * unprim_frac; 

            if unprim_onestate
                SS_empty = [steadyvesicles_reshaped(1) steadyvesicles_reshaped(num_states+1) steadyvesicles_reshaped(2*num_states+1)] * (Caunprim_rate/Caprim_rate);
            end

        elseif SS_PM == 1

            steady_slowstates = zeros(3*(m_max+1),1);

            for yy = 1:3 %Activation states
                for zz = 1:(m_max+1)
                    steady_slowstates((yy-1)*(m_max+1) + zz) = sum(steadyvesicles_reshaped((((yy-1)*num_states)+zz):(m_max+1):((yy-1)*num_states + (n_max*(m_max+1)+zz))));
                end
            end

            SS_empty = steady_slowstates * unprim_frac;

        end            

        if ~unprim_onestate
            steadyvesicles_reshaped(1:num_states_total) = steadyvesicles_reshaped(1:num_states_total) * prim_frac;
        end

        steadyvesicles_reshaped((num_states_total+1):(num_states_total + num_empty_fus -1)) = SS_empty;  
        steadyvesicles_reshaped = steadyvesicles_reshaped* (num_ves_in_steady/sum(steadyvesicles_reshaped));

    elseif u_val ~= 1

        steady_priming = calculate_steady_state_priming(par, Ca_bas); %steady_primin = [P, R0, R1, R2]; Total sum: 1.

            par_alt = par;
            par_alt([6:7 14]) = 0;
            steady_reshaped_first_sensor = calculate_steady_state(par_alt, Ca_bas, num_ves_in_steady); %First sensor steady


        for kk = 1:(m_max+1)
            inds = kk:(m_max+1):(((n_max+1)-1)*(m_max+1)+kk);
            steadyvesicles_reshaped(inds) = steady_reshaped_first_sensor'*steady_priming(kk+1);
        end
            steadyvesicles_reshaped(num_states_total+1) = num_ves_in_steady*steady_priming(1);
    end

end
    

time_cells = cell(1,num_bins);
vesicles_cells = cell(1,num_bins);
ves_avail = zeros(num_bins, length(time_vector));
    
if pVr2_hack == 0
    parfor qq = 1:num_bins
        %disp(['BEFORE simulation ' num2str(qq)])
        Ca_R_sim = Ca_R_vesicles(:,qq);
        dist = dists_vesicles(qq);
        [t, vesicles] = ode15s(@(t, input_pop)dual_sensor_ODE_FUSION(act_model_type, SS_PM, Ca_prim_type, Ca_R_sim, dist, input_pop, par, Ca_time_vesicles, t), [0 simulation_time], steadyvesicles_reshaped);
        %disp(['AFTER simulation ' num2str(qq)])
        time_cells{qq} = t;
        vesicles_cells{qq} = vesicles;
    end
elseif pVr2_hack == 1
    parfor qq = 1:num_bins
        %disp(['BEFORE simulation ' num2str(qq)])
        simulation_time1 = 0.01; %First round: 10 ms
        Ca_R_sim = Ca_R_vesicles(:,qq);
        dist = dists_vesicles(qq);
        [t1, vesicles1] = ode15s(@(t, input_pop)dual_sensor_ODE_FUSION(act_model_type, SS_PM, Ca_prim_type, Ca_R_sim, dist, input_pop, par, Ca_time_vesicles, t), [0 simulation_time1], steadyvesicles_reshaped);
        %disp(['AFTER simulation ' num2str(qq)])
        
        t2_init = t1(end)+1e-6; 
        
        steadyvesicles_reshaped2 = vesicles1(end,:);
        steadyvesicles_reshaped2(1:num_states_total) = steadyvesicles_reshaped2(1:num_states_total)*(num_ves_in_steady/sum(steadyvesicles_reshaped2(1:num_states_total)));
        steadyvesicles_reshaped2((num_states_total+1):end) = 0;
        [t2, vesicles2] = ode15s(@(t, input_pop)dual_sensor_ODE_FUSION(act_model_type, SS_PM, Ca_prim_type, Ca_R_sim, dist, input_pop, par, Ca_time_vesicles, t), [t2_init simulation_time], steadyvesicles_reshaped2);
        
        t = [t1; t2];
        vesicles = [vesicles1; vesicles2];
        
        
        time_cells{qq} = t;
        vesicles_cells{qq} = vesicles;
    end
end


for j = 1:num_bins
    time_outcome = time_cells{j};
    ves_outcome = vesicles_cells{j};
    ves_interp = interp1(time_outcome, ves_outcome, time_vector);
    allstatesinallbins_cell{j} = ves_interp;

    states_all = states_all + ves_interp;
    ves_avail(j,:) = sum(ves_interp(:,1:num_states),2)/num_ves_in_steady;

    fused_states(:,j) = fused_states(:,j) + ves_interp(:,fuse_state);
    empty_states_1(:,j) = ves_interp(:,3*(n_max+1)*(m_max+1)+1);
    empty_states_2(:,j) = ves_interp(:,3*(n_max+1)*(m_max+1)+2);
    empty_states_3(:,j) = ves_interp(:,3*(n_max+1)*(m_max+1)+3);
    fused_states_EPSC = fused_states_EPSC+squeeze(fused_states(:,j));
end
[fast_states, slow_states, total_fast_states, total_slow_states] = determine_total_slowfast_DET(time_vector, states_all, par);
fused_states_EPSC = fused_states_EPSC';






if pVr2_hack == 0
    [~, mEPSC_current] = smooth_mEPSC;
    mEPSC_init = mEPSC_current;
    mEPSC = -(mEPSC_init/min(mEPSC_init))*size_of_mini;


    fused_ves_vec = [0 diff(fused_states_EPSC)];
    EPSC = conv(fused_ves_vec, mEPSC);
    EPSC = EPSC(1:length(time_vector));
else
    fused_ves_vec = NaN;
    EPSC = NaN;
    EPSC = NaN;
end

pVr1_bins_10ms = squeeze(fused_states(10000,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr1_bins_3_5ms = squeeze(fused_states(3500,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr1_bins_4ms = squeeze(fused_states(4000,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr1_bins_4_5ms = squeeze(fused_states(4500,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr1_bins_5ms = squeeze(fused_states(5000,:))/(sum(steadyvesicles_reshaped(1:num_states)));

pVr2_bins_20ms = squeeze(fused_states(20000,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr2_bins_13_5ms = squeeze(fused_states(13500,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr2_bins_14ms = squeeze(fused_states(14000,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr2_bins_14_5ms = squeeze(fused_states(14500,:))/(sum(steadyvesicles_reshaped(1:num_states)));
pVr2_bins_15ms = squeeze(fused_states(15000,:))/(sum(steadyvesicles_reshaped(1:num_states)));






pVr1 = fused_states_EPSC(10000)/num_ves;


if pVr2_hack == 0
    [peak1, peak2, peaksum, ppr, peak2_extr, ppr_extr, ttp2080, peaksum_extr] = determine_peak_ppr(time_vector', EPSC', stim_freq);
else
    [time_vector, EPSC, peak1, peak2, ppr, states_all, pVr1, ttp2080,  peak2_extr, ppr_extr, fast_states, slow_states, total_fast_states, total_slow_states, peaksum_extr] = deal(NaN);
end
