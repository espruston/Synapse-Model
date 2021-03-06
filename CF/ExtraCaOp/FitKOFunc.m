function err = FitKOFunc(x)

CFData3KO = load('..\CFData3KO.mat').CFData3KO;
CF1HzCa = load('..\CF1HzCa.mat').Ca_sim;
CF10HzCa = load('..\CF10HzCa.mat').Ca_sim;
CF20HzCa = load('..\CF20HzCa.mat').Ca_sim;
CF50HzCa = load('..\CF50HzCa.mat').Ca_sim;

p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);

C_Ca = x(5);

Ca_rest = 50e-9;

%Define train stims
stimulus_times_1 = linspace(0,1000*19,20); %100 stims 1hz
max_time_1 = stimulus_times_1(end) + stimulus_times_1(2)*3;

stimulus_times_10 = linspace(0,100*19,20); %100 stims 10hz
max_time_10 = stimulus_times_10(end) + stimulus_times_10(2)*3;

stimulus_times_20 = linspace(0,50*19,20);
max_time_20 = stimulus_times_20(end) + stimulus_times_20(2)*3;

stimulus_times_50 = linspace(0,20*19,20);
max_time_50 = stimulus_times_50(end) + stimulus_times_50(2)*3;

t_SS = 10000; %ms
state_0 = [1; 0; reserve_size]; %[empty pool; bound pool; reserve pool]

%Syt3 KO
[t0,state] = ode15s(@(t,state) dSS(t, state, k_docking, k_undocking, Ca_rest*C_Ca, reserve_size), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, CF1HzCa*C_Ca, SS, reserve_size);
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, CF10HzCa*C_Ca, SS, reserve_size);
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, CF20HzCa*C_Ca, SS, reserve_size);
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, CF50HzCa*C_Ca, SS, reserve_size);
hz_50 = Fused_50/Fused_50(1);

err_3KO = sqrt(sum((hz_1 - CFData3KO(1:20,1)).^2 + (hz_10 - CFData3KO(1:20,2)).^2 + (hz_20 - CFData3KO(1:20,3)).^2 + (hz_50(1:20) - CFData3KO(1:20,4)).^2) + sum(10*(hz_50(21:30) - CFData3KO(21:30,4)).^2));

err = err_3KO;

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, Ca, SS, reserve_size)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release, -pre_stim(2)*p_release, 0];
        Fused(i) = pre_stim(2)*p_release;  
        [t,out] = ode15s(@(t,state) dState(t, state, k_docking, k_undocking, Ca(1,:), reserve_size), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [t,out] = ode15s(@(t,state) dState(t, state, k_docking, k_undocking, Ca(i,:), reserve_size), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,Ca,reserve_size)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+Ca) + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1);

end

function dydt = dState(t,state,k_docking,k_undocking,Ca,reserve_size)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+Ca(round(t/.01)+1)) + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1);

end

end