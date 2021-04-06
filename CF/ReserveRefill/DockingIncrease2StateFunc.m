function cost = DockingIncrease2StateFunc(x)
disp(['Testing x = [', num2str(x), ']'])
CFDataWT = load('..\CFDataWT.mat').CFDataWT;
CFDataKO = load('..\CFData3KO.mat').CFData3KO;
CF1HzSyt3 = load('..\CF1HzSyt3.mat').Syt3;
CF10HzSyt3 = load('..\CF10HzSyt3.mat').Syt3;
CF20HzSyt3 = load('..\CF20HzSyt3.mat').Syt3;
CF50HzSyt3 = load('..\CF50HzSyt3.mat').Syt3;

p_release_dock = x(1); 
p_release_tether = x(2);
k_docking = x(3); 
k_undocking = x(4); 
k_tether = x(5);
k_untether = x(6);
reserve_size = x(7);
k_refill = x(8);

C_3 = x(9); %multiplicative factor of syt3 effect

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
state_0 = [1; 0; 0; reserve_size]; %[empty pool; docked pool; tethered pool; reserve pool]

%WT
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, C_3*CF1HzSyt3(1)), [0 t_SS], state_0);

SS_WT = state(end,:);

[~, ~, Fused_1_WT] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*CF1HzSyt3);
hz_1_WT = Fused_1_WT/Fused_1_WT(1);

[~, ~, Fused_10_WT] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*CF10HzSyt3);
hz_10_WT = Fused_10_WT/Fused_10_WT(1);

[~, ~, Fused_20_WT] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*CF20HzSyt3);
hz_20_WT = Fused_20_WT/Fused_20_WT(1);

[~, ~, Fused_50_WT] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*CF50HzSyt3);
hz_50_WT = Fused_50_WT/Fused_50_WT(1);

err_WT = sqrt(2*sum((hz_1_WT(1:5) - CFDataWT(1:5,1)).^2 + (hz_10_WT(1:5) - CFDataWT(1:5,2)).^2 + (hz_20_WT(1:5) - CFDataWT(1:5,3)).^2 + (hz_50_WT(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_1_WT(6:20) - CFDataWT(6:20,1)).^2 + (hz_10_WT(6:20) - CFDataWT(6:20,2)).^2 + (hz_20_WT(6:20) - CFDataWT(6:20,3)).^2 + (hz_50_WT(6:20) - CFDataWT(6:20,4)).^2) + 10*sum((hz_50_WT(21:30) - CFDataWT(21:30,4)).^2));

%KO
C_3 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, C_3*CF1HzSyt3(1)), [0 t_SS], state_0);

SS_KO = state(end,:);

[~, ~, Fused_1_KO] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_KO, reserve_size, k_refill, C_3*CF1HzSyt3);
hz_1_KO = Fused_1_KO/Fused_1_KO(1);

[~, ~, Fused_10_KO] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_KO, reserve_size, k_refill, C_3*CF10HzSyt3);
hz_10_KO = Fused_10_KO/Fused_10_KO(1);

[~, ~, Fused_20_KO] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_KO, reserve_size, k_refill, C_3*CF20HzSyt3);
hz_20_KO = Fused_20_KO/Fused_20_KO(1);

[~, ~, Fused_50_KO] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_KO, reserve_size, k_refill, C_3*CF50HzSyt3);
hz_50_KO = Fused_50_KO/Fused_50_KO(1);

err_KO = sqrt(2*sum((hz_1_KO(1:5) - CFDataKO(1:5,1)).^2 + (hz_10_KO(1:5) - CFDataKO(1:5,2)).^2 + (hz_20_KO(1:5) - CFDataKO(1:5,3)).^2 + (hz_50_KO(1:5) - CFDataKO(1:5,4)).^2) + sum((hz_1_KO(6:20) - CFDataKO(6:20,1)).^2 + (hz_10_KO(6:20) - CFDataKO(6:20,2)).^2 + (hz_20_KO(6:20) - CFDataKO(6:20,3)).^2 + (hz_50_KO(6:20) - CFDataKO(6:20,4)).^2) + 10*sum((hz_50_KO(21:30) - CFDataKO(21:30,4)).^2));

err = (err_WT + err_KO)/2;

cost = err + abs(err_WT - err_KO)/10;

disp(['Cost = ', num2str(cost), ', average error = ', num2str(err), ', WT error = ', num2str(err_WT), ', KO error = ', num2str(err_KO)])


function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS, reserve_size, k_refill, Syt3)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release_dock + pre_stim(3)*p_release_tether, -pre_stim(2)*p_release_dock, -pre_stim(3)*p_release_tether, 0];
        Fused(i) = pre_stim(2)*p_release_dock + pre_stim(3)*p_release_tether;  
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, Syt3(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [~,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, Syt3(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release_dock + pre_stim(3)*p_release_tether;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,k_tether,k_untether,reserve_size,k_refill,Syt3)

dydt(1,1) = -state(1)*(state(4)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1) - state(2)*k_tether*(1+Syt3) + state(3)*k_untether;
dydt(3,1) = state(2)*k_tether*(1+Syt3) - state(3)*k_untether;
dydt(4,1) = dydt(1,1) + (reserve_size-state(4))*k_refill;

end

function dydt = dState(t,state,k_docking,k_undocking,k_tether,k_untether,reserve_size,k_refill,Syt3)

dydt(1,1) = -state(1)*(state(4)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1) - state(2)*k_tether*(1+Syt3(round(t/.1)+1)) + state(3)*k_untether;
dydt(3,1) = state(2)*k_tether*(1+Syt3(round(t/.1)+1)) - state(3)*k_untether;
dydt(4,1) = dydt(1,1) + (reserve_size-state(4))*k_refill;

end
end