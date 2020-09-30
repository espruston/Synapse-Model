function err = MaturationCDRFunc(x)
%MATURATIONCDRFUNC Summary of this function goes here
%   Detailed explanation goes here

p_mature = x(1); 
p_immature = x(2); 
k_docking = x(3); 
k_undocking = x(4); 
k_maturation = x(5); 
k_dematuration = x(6); 

data = matfile('WT_data.mat').WT_data;
hz_1_data = data(:,1);
hz_10_data = data(:,2);
hz_20_data = data(:,3);
hz_50_data = data(:,4);
hz_100_data = data(:,5);
hz_200_data = data(:,6);

k_on_3 = 3e5; %M^-1ms^-1 Hui
k_off_3 = 0.05; %ms^-1  Hui
k_on_7 = 7.333e3; %Knight
k_off_7 = 1.1e-2;
C_3 = 1;
Ca_rest = 5e-8; %M

t_SS = 10000; %ms
%ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [1; 0; 0; 0; 0]; %[empty; immature; mature; syt3; syt7]

[t0,state] = ode15s(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_rest), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, hz_1, ~, ~, hz_10, ~, ~, hz_20, ~, ~, hz_50, ~, ~, hz_100, ~, ~, hz_200] = test6(p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);

err = sqrt(sum((hz_1 - hz_1_data).^2 + (hz_10 - hz_10_data).^2 + (hz_20 - hz_20_data).^2 + (hz_50 - hz_50_data).^2 + (hz_100 - hz_100_data).^2 + (hz_200 - hz_200_data).^2));
  
function [ts, state, Fused_im, Fused_m, Ca_sim] = stim_sim(stimulus_times, max_time, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS)

    k_on_3 = 3e5; %M^-1ms^-1 Hui
    k_off_3 = 0.05; %ms^-1  Hui
    k_on_7 = 7.333e3; %Knight
    k_off_7 = 1.1e-2;
    C_3 = 1;
    C_7 = 0;
    delta_t = 1e-2; %ms
    
    state = SS;
    Ca_sim = create_Ca_signal(stimulus_times, max_time);
    %[syt3, syt7] = syt_sim(Ca_sim);
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused_im = zeros(length(stimulus_times),1);
    Fused_m = zeros(length(stimulus_times),1);

    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_immature+pre_stim(3)*p_mature, -pre_stim(2)*p_immature, -pre_stim(3)*p_mature, 0, 0];
        Fused_im(i) = pre_stim(2)*p_immature;
        Fused_m(i) = pre_stim(3)*p_mature;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_sim), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end

    state = [SS; state];
    ts = [delta_t; ts];
end


function Ca_sim = create_Ca_signal(stimulus_times, max_time)
    
    FWHM = .34; %Local calcium full width half maximum ms
    sigma = FWHM/2.35; %variance
    mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

    Ca_rest = 5e-8; %M
    Ca_spike = 2e-5; %M
    Ca_residual = 250e-9; %M
    T_Ca_decay = 40; %ms
    
    delta_t = 1e-2;
    
    
    ts = linspace(0,max_time,max_time/delta_t + 1);
    Ca_sim = zeros(1,length(ts));
    Ca_sim = Ca_sim + Ca_rest;

    for t = 1:length(stimulus_times) %simulate calcium influx

        spike_start_index = round(stimulus_times(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times(t)+mu)/delta_t) + 1;

        Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one

        Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

    end

end


function [Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50, Fused_im_100, Fused_m_100, hz_100, Fused_im_200, Fused_m_200, hz_200] = test6(p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS)
    
    stimulus_times_1 = linspace(0,1000*99,100); %100 stims 1hz
    stimulus_times_10 = linspace(0,100*99,100); %100 stims 10hz
    stimulus_times_20 = linspace(0,50*99,100); 
    stimulus_times_50 = linspace(0,20*99,100);
    stimulus_times_100 = linspace(0,10*99,100);
    stimulus_times_200 = linspace(0,5*99,100);
    

    max_time_1 = stimulus_times_1(end) + stimulus_times_1(2)*30;
    [~, ~, Fused_im_1, Fused_m_1,~] = stim_sim(stimulus_times_1, max_time_1, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);
    Fused_1 = Fused_im_1 + Fused_m_1;
    hz_1 = Fused_1/Fused_1(1);
    
    max_time_10 = stimulus_times_10(end) + stimulus_times_10(2)*30;
    [~, ~, Fused_im_10, Fused_m_10,~] = stim_sim(stimulus_times_10, max_time_10, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);
    Fused_10 = Fused_im_10 + Fused_m_10;
    hz_10 = Fused_10/Fused_10(1);
    
    max_time_20 = stimulus_times_20(end) + stimulus_times_20(2)*30;
    [~, ~, Fused_im_20, Fused_m_20,~] = stim_sim(stimulus_times_20, max_time_20, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);
    Fused_20 = Fused_im_20 + Fused_m_20;
    hz_20 = Fused_20/Fused_20(1);

    max_time_50 = stimulus_times_50(end) + stimulus_times_50(2)*30;
    [~, ~, Fused_im_50, Fused_m_50, ~] = stim_sim(stimulus_times_50, max_time_50, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);
    Fused_50 = Fused_im_50 + Fused_m_50;
    hz_50 = Fused_50/Fused_50(1);
    
    max_time_100 = stimulus_times_100(end) + stimulus_times_100(2)*100;
    [~, ~, Fused_im_100, Fused_m_100, ~] = stim_sim(stimulus_times_100, max_time_100, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);
    Fused_100 = Fused_im_100 + Fused_m_100;
    hz_100 = Fused_100/Fused_100(1);
    
    max_time_200 = stimulus_times_200(end) + stimulus_times_200(2)*30;
    [~, ~, Fused_im_200, Fused_m_200, ~] = stim_sim(stimulus_times_200, max_time_200, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, SS);
    Fused_200 = Fused_im_200 + Fused_m_200;
    hz_200 = Fused_200/Fused_200(1);

end

function dydt = dSS(~,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_rest)
    
    k_docking = k_docking*(1+state(4)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    dydt(4,1) = (1-state(4))*k_on_3*Ca_rest - state(4)*k_off_3;
    dydt(5,1) = (1-state(5))*k_on_7*Ca_rest - state(5)*k_off_7;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_sim)
    
    Ca = Ca_sim(round(t/.01)+1);
    k_docking = k_docking*(1+state(4)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    dydt(4,1) = (1-state(4))*k_on_3*Ca - state(4)*k_off_3;
    dydt(5,1) = (1-state(5))*k_on_7*Ca - state(5)*k_off_7;
    
end

end

