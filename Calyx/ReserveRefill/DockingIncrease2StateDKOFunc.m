function cost = DockingIncrease2StateDKOFunc(x)
disp(['Testing x = [', num2str(x), ']'])

CalyxDataDKO = load('..\CalyxDataDKO.mat').CalyxDataDKO;

%x = [0.01,0.25,0.00605823,0.0731115,0.0250303,0.00192173,11.5335,9.39172e-05]; %

p_release_dock = x(1); 
p_release_tether = x(2);
k_docking = x(3); 
k_undocking = x(4);
k_tether = x(5);
k_untether = x(6);
reserve_size = x(7);
k_refill = x(8);

%Define train stims
stimulus_times_1 = linspace(0,1000*99,100); %100 stims 1hz
max_time_1 = stimulus_times_1(end) + stimulus_times_1(2)*3;

stimulus_times_10 = linspace(0,100*99,100); %100 stims 10hz
max_time_10 = stimulus_times_10(end) + stimulus_times_10(2)*3;

stimulus_times_20 = linspace(0,50*99,100);
max_time_20 = stimulus_times_20(end) + stimulus_times_20(2)*3;

stimulus_times_50 = linspace(0,20*99,100);
max_time_50 = stimulus_times_50(end) + stimulus_times_50(2)*3;

stimulus_times_100 = linspace(0,10*99,100); %100 stims 10hz
max_time_100 = stimulus_times_100(end) + stimulus_times_100(2)*3;

stimulus_times_200 = linspace(0,5*99,100);
max_time_200 = stimulus_times_200(end) + stimulus_times_200(2)*3;

t_SS = 10000; %ms
state_0 = [1; 0; 0; reserve_size]; %[empty pool; docked pool; tethered pool, reserve pool]

%DKO
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill), [0 t_SS], state_0);

SS_DKO = state(end,:);

[~, ~, Fused_1_DKO] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_1_DKO = Fused_1_DKO/Fused_1_DKO(1);

[~, ~, Fused_10_DKO] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_10_DKO = Fused_10_DKO/Fused_10_DKO(1);

[~, ~, Fused_20_DKO] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_20_DKO = Fused_20_DKO/Fused_20_DKO(1);

[~, ~, Fused_50_DKO] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_50_DKO = Fused_50_DKO/Fused_50_DKO(1);

[~, ~, Fused_100_DKO] = stim_sim(stimulus_times_100, max_time_100, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_100_DKO = Fused_100_DKO/Fused_100_DKO(1);

[~, ~, Fused_200_DKO] = stim_sim(stimulus_times_200, max_time_200, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_200_DKO = Fused_200_DKO/Fused_200_DKO(1);

err_DKO = sqrt(2*sum((hz_1_DKO(1:10) - CalyxDataDKO(1:10,1)).^2 + (hz_10_DKO(1:10) - CalyxDataDKO(1:10,2)).^2 + (hz_20_DKO(1:10) - CalyxDataDKO(1:10,3)).^2 + (hz_50_DKO(1:10) - CalyxDataDKO(1:10,4)).^2 + (hz_100_DKO(1:10) - CalyxDataDKO(1:10,5)).^2 + (hz_200_DKO(1:10) - CalyxDataDKO(1:10,6)).^2) + sum((hz_1_DKO(11:100) - CalyxDataDKO(11:100,1)).^2 + (hz_10_DKO(11:100) - CalyxDataDKO(11:100,2)).^2 + (hz_20_DKO(11:100) - CalyxDataDKO(11:100,3)).^2 + (hz_50_DKO(11:100) - CalyxDataDKO(11:100,4)).^2 + (hz_100_DKO(11:100) - CalyxDataDKO(11:100,5)).^2 + (hz_200_DKO(11:100) - CalyxDataDKO(11:100,6)).^2) + 10*sum((hz_10_DKO(101:111) - CalyxDataDKO(101:111,2)).^2 + (hz_20_DKO(101:111) - CalyxDataDKO(101:111,3)).^2 + (hz_50_DKO(101:111) - CalyxDataDKO(101:111,4)).^2 + (hz_100_DKO(101:111) - CalyxDataDKO(101:111,5)).^2 + (hz_200_DKO(101:111) - CalyxDataDKO(101:111,6)).^2));

err = err_DKO;

cost = err + abs(err - err_DKO)/5;

disp(['Cost = ', num2str(cost), ', DKO error = ', num2str(err_DKO)])

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS, reserve_size, k_refill)

    
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay stim_delay(1)];
    if ~isequal(stimulus_times, linspace(0,1000*99,100)) %1hz
        Recovery = [10 20 50 100 200 500 1000 2000 5000 10000 20000];
        Recovery = [10 diff(Recovery)];
        stim_delay = [stim_delay Recovery];
    end

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release_dock + pre_stim(3)*p_release_tether, -pre_stim(2)*p_release_dock, -pre_stim(3)*p_release_tether, 0];
        Fused(i) = pre_stim(2)*p_release_dock + pre_stim(3)*p_release_tether;
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
end

function dydt = dSS(~,state,k_docking,k_undocking,k_tether,k_untether,reserve_size,k_refill)

dydt(1,1) = -state(1)*(state(4)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1) - state(2)*k_tether + state(3)*k_untether;
dydt(3,1) = state(2)*k_tether - state(3)*k_untether;
dydt(4,1) = dydt(1,1) + (reserve_size-state(4))*k_refill;

end

function dydt = dState(~,state,k_docking,k_undocking,k_tether,k_untether,reserve_size,k_refill)

dydt(1,1) = -state(1)*(state(4)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1) - state(2)*k_tether + state(3)*k_untether;
dydt(3,1) = state(2)*k_tether - state(3)*k_untether;
dydt(4,1) = dydt(1,1) + (reserve_size-state(4))*k_refill;

end

end