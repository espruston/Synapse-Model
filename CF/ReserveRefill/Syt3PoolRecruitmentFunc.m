function cost = Syt3PoolRecruitmentFunc(x)
disp(['Testing x = [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), ',', num2str(x(7)), ']'])

CFDataWT = load('..\CFDataWT.mat').CFDataWT;
CFDataKO = load('..\CFData3KO.mat').CFData3KO;
CF1HzSyt3 = load('..\CF1HzSyt3.mat').Syt3;
CF10HzSyt3 = load('..\CF10HzSyt3.mat').Syt3;
CF20HzSyt3 = load('..\CF20HzSyt3.mat').Syt3;
CF50HzSyt3 = load('..\CF50HzSyt3.mat').Syt3;

p_release = x(1); 
k_docking = x(2);
k_undocking = x(3);
reserve_size = x(4);
k_refill = x(5);

Syt_pool_size = x(6);

C_3 = x(7);

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
state_0 = [1-Syt_pool_size; 0; Syt_pool_size; 0; reserve_size]; %[empty pool; bound pool; empty syt pool; bound syt pool; reserve pool]

%WT
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, C_3*CF1HzSyt3(1)), [0 t_SS], state_0);

SS_WT = state(end,:);

[~, ~, Fused_1_WT] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*CF1HzSyt3);
hz_1_WT = Fused_1_WT/Fused_1_WT(1);

[~, ~, Fused_10_WT] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*CF10HzSyt3);
hz_10_WT = Fused_10_WT/Fused_10_WT(1);

[~, ~, Fused_20_WT] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*CF20HzSyt3);
hz_20_WT = Fused_20_WT/Fused_20_WT(1);

[~, ~, Fused_50_WT] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*CF50HzSyt3);
hz_50_WT = Fused_50_WT/Fused_50_WT(1);

err_WT = sqrt(2*sum((hz_1_WT(1:5) - CFDataWT(1:5,1)).^2 + (hz_10_WT(1:5) - CFDataWT(1:5,2)).^2 + (hz_20_WT(1:5) - CFDataWT(1:5,3)).^2 + (hz_50_WT(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_1_WT(6:20) - CFDataWT(6:20,1)).^2 + (hz_10_WT(6:20) - CFDataWT(6:20,2)).^2 + (hz_20_WT(6:20) - CFDataWT(6:20,3)).^2 + (hz_50_WT(6:20) - CFDataWT(6:20,4)).^2) + 10*sum((hz_50_WT(21:30) - CFDataWT(21:30,4)).^2));

%KO
C_3 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, C_3*CF1HzSyt3(1)), [0 t_SS], state_0);

SS_KO = state(end,:);

[~, ~, Fused_1_KO] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, Syt_pool_size, C_3*CF1HzSyt3);
hz_1_KO = Fused_1_KO/Fused_1_KO(1);

[~, ~, Fused_10_KO] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, Syt_pool_size, C_3*CF10HzSyt3);
hz_10_KO = Fused_10_KO/Fused_10_KO(1);

[~, ~, Fused_20_KO] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, Syt_pool_size, C_3*CF20HzSyt3);
hz_20_KO = Fused_20_KO/Fused_20_KO(1);

[~, ~, Fused_50_KO] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, Syt_pool_size, C_3*CF50HzSyt3);
hz_50_KO = Fused_50_KO/Fused_50_KO(1);

err_KO = sqrt(2*sum((hz_1_KO(1:5) - CFDataKO(1:5,1)).^2 + (hz_10_KO(1:5) - CFDataKO(1:5,2)).^2 + (hz_20_KO(1:5) - CFDataKO(1:5,3)).^2 + (hz_50_KO(1:5) - CFDataKO(1:5,4)).^2) + sum((hz_1_KO(6:20) - CFDataKO(6:20,1)).^2 + (hz_10_KO(6:20) - CFDataKO(6:20,2)).^2 + (hz_20_KO(6:20) - CFDataKO(6:20,3)).^2 + (hz_50_KO(6:20) - CFDataKO(6:20,4)).^2) + 10*sum((hz_50_KO(21:30) - CFDataKO(21:30,4)).^2));

err = (err_WT + err_KO)/2;

cost = err + abs(err_WT - err_KO)/10;

disp(['Cost = ', num2str(cost), ', average error = ', num2str(err), ', WT error = ', num2str(err_WT), ', KO error = ', num2str(err_KO)])

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, SS, reserve_size, k_refill, Syt_pool_size, Syt3)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release, -pre_stim(2)*p_release, pre_stim(4)*p_release, -pre_stim(4)*p_release, 0];
        Fused(i) = pre_stim(2)*p_release;  
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, Syt3(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [~,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, Syt3(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release + pre_stim(4)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,reserve_size,k_refill,Syt_pool_size,Syt3)

dydt(1,1) = -(state(1)/(1-Syt_pool_size))*(state(5)/reserve_size)*k_docking + (state(2)/(1-Syt_pool_size))*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = -(state(3)/Syt_pool_size)*(state(5)/reserve_size)*k_docking*Syt3 + (state(4)/Syt_pool_size)*k_undocking;
dydt(4,1) = -dydt(3,1);
dydt(5,1) = dydt(1,1) + dydt(3,1) + (reserve_size-state(5))*k_refill;
%dydt(6,1) = (1-state(6))*k_on_3*(Ca_rest^3.6/(Ca_rest^3.6+7e-6^3.6)) - state(6)*k_off_3;
%dydt(6,1) = (1-state(6))*k_on_3*Ca_rest - state(6)*k_off_3;

end

function dydt = dState(t,state,k_docking,k_undocking,reserve_size,k_refill,Syt_pool_size,Syt3)

dydt(1,1) = -(state(1)/(1-Syt_pool_size))*(state(5)/reserve_size)*k_docking + (state(2)/(1-Syt_pool_size))*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = -(state(3)/Syt_pool_size)*(state(5)/reserve_size)*k_docking*Syt3(round(t/.1)+1) + (state(4)/Syt_pool_size)*k_undocking;
dydt(4,1) = -dydt(3,1);
dydt(5,1) = dydt(1,1) + dydt(3,1) + (reserve_size-state(5))*k_refill;
%dydt(6,1) = (1-state(6))*k_on_3*(Ca(round(t/.1)+1)^3.6/(Ca(round(t/.1)+1)^3.6+7e-6^3.6)) - state(6)*k_off_3;
%dydt(6,1) = (1-state(6))*k_on_3*(Ca(round(t/.1)+1)) - state(6)*k_off_3;

end

end