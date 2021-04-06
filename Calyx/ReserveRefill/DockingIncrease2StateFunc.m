function cost = DockingIncrease2StateFunc(x)
disp(['Testing x = [', num2str(x), ']'])
CalyxDataWT = load('..\CalyxDataWT.mat').CalyxDataWT;
CalyxData3KO = load('..\CalyxData3KO.mat').CalyxData3KO;
CalyxData7KO = load('..\CalyxData7KO.mat').CalyxData7KO;
CalyxDataDKO = load('..\CalyxDataDKO.mat').CalyxDataDKO;
Calyx1HzSyt3 = load('..\Calyx1HzSyt3.mat').Syt3;
Calyx10HzSyt3 = load('..\Calyx10HzSyt3.mat').Syt3;
Calyx20HzSyt3 = load('..\Calyx20HzSyt3.mat').Syt3;
Calyx50HzSyt3 = load('..\Calyx50HzSyt3.mat').Syt3;
Calyx100HzSyt3 = load('..\Calyx100HzSyt3.mat').Syt3;
Calyx200HzSyt3 = load('..\Calyx200HzSyt3.mat').Syt3;
Calyx1HzSyt7 = load('..\Calyx1HzSyt7.mat').Syt7;
Calyx10HzSyt7 = load('..\Calyx10HzSyt7.mat').Syt7;
Calyx20HzSyt7 = load('..\Calyx20HzSyt7.mat').Syt7;
Calyx50HzSyt7 = load('..\Calyx50HzSyt7.mat').Syt7;
Calyx100HzSyt7 = load('..\Calyx100HzSyt7.mat').Syt7;
Calyx200HzSyt7 = load('..\Calyx200HzSyt7.mat').Syt7;

p_release_dock = x(1); 
p_release_tether = x(2);
k_docking = x(3); 
k_undocking = x(4);
k_tether = x(5);
k_untether = x(6);
reserve_size = x(7);
k_refill = x(8);

C_3 = x(9); %additive factor of syt3 effect
C_7 = x(10); %syt7 effect

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
state_0 = [1; 0; reserve_size]; %[empty pool; bound pool; reserve pool]

%WT
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_WT = state(end,:);

[~, ~, Fused_1_WT] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_WT = Fused_1_WT/Fused_1_WT(1);

[~, ~, Fused_10_WT] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_WT = Fused_10_WT/Fused_10_WT(1);

[~, ~, Fused_20_WT] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_WT = Fused_20_WT/Fused_20_WT(1);

[~, ~, Fused_50_WT] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_WT = Fused_50_WT/Fused_50_WT(1);

[~, ~, Fused_100_WT] = stim_sim(stimulus_times_100, max_time_100, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_WT = Fused_100_WT/Fused_100_WT(1);

[~, ~, Fused_200_WT] = stim_sim(stimulus_times_200, max_time_200, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_WT, reserve_size, k_refill, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_WT = Fused_200_WT/Fused_200_WT(1);

err_WT = sqrt(2*sum((hz_1_WT(1:10) - CalyxDataWT(1:10,1)).^2 + (hz_10_WT(1:10) - CalyxDataWT(1:10,2)).^2 + (hz_20_WT(1:10) - CalyxDataWT(1:10,3)).^2 + (hz_50_WT(1:10) - CalyxDataWT(1:10,4)).^2 + (hz_100_WT(1:10) - CalyxDataWT(1:10,5)).^2 + (hz_200_WT(1:10) - CalyxDataWT(1:10,6)).^2) + sum((hz_1_WT(11:100) - CalyxDataWT(11:100,1)).^2 + (hz_10_WT(11:100) - CalyxDataWT(11:100,2)).^2 + (hz_20_WT(11:100) - CalyxDataWT(11:100,3)).^2 + (hz_50_WT(11:100) - CalyxDataWT(11:100,4)).^2 + (hz_100_WT(11:100) - CalyxDataWT(11:100,5)).^2 + (hz_200_WT(11:100) - CalyxDataWT(11:100,6)).^2) + 10*sum((hz_10_WT(101:111) - CalyxDataWT(101:111,2)).^2 + (hz_20_WT(101:111) - CalyxDataWT(101:111,3)).^2 + (hz_50_WT(101:111) - CalyxDataWT(101:111,4)).^2 + (hz_100_WT(101:111) - CalyxDataWT(101:111,5)).^2 + (hz_200_WT(101:111) - CalyxDataWT(101:111,6)).^2));

%3KO
C_3 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_3KO = state(end,:);

[~, ~, Fused_1_3KO] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_3KO, reserve_size, k_refill, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_3KO = Fused_1_3KO/Fused_1_3KO(1);

[~, ~, Fused_10_3KO] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_3KO, reserve_size, k_refill, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_3KO = Fused_10_3KO/Fused_10_3KO(1);

[~, ~, Fused_20_3KO] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_3KO, reserve_size, k_refill, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_3KO = Fused_20_3KO/Fused_20_3KO(1);

[~, ~, Fused_50_3KO] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_3KO, reserve_size, k_refill, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_3KO = Fused_50_3KO/Fused_50_3KO(1);

[~, ~, Fused_100_3KO] = stim_sim(stimulus_times_100, max_time_100, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_3KO, reserve_size, k_refill, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_3KO = Fused_100_3KO/Fused_100_3KO(1);

[~, ~, Fused_200_3KO] = stim_sim(stimulus_times_200, max_time_200, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_3KO, reserve_size, k_refill, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_3KO = Fused_200_3KO/Fused_200_3KO(1);

err_3KO = sqrt(2*sum((hz_1_3KO(1:10) - CalyxData3KO(1:10,1)).^2 + (hz_10_3KO(1:10) - CalyxData3KO(1:10,2)).^2 + (hz_20_3KO(1:10) - CalyxData3KO(1:10,3)).^2 + (hz_50_3KO(1:10) - CalyxData3KO(1:10,4)).^2 + (hz_100_3KO(1:10) - CalyxData3KO(1:10,5)).^2 + (hz_200_3KO(1:10) - CalyxData3KO(1:10,6)).^2) + sum((hz_1_3KO(11:100) - CalyxData3KO(11:100,1)).^2 + (hz_10_3KO(11:100) - CalyxData3KO(11:100,2)).^2 + (hz_20_3KO(11:100) - CalyxData3KO(11:100,3)).^2 + (hz_50_3KO(11:100) - CalyxData3KO(11:100,4)).^2 + (hz_100_3KO(11:100) - CalyxData3KO(11:100,5)).^2 + (hz_200_3KO(11:100) - CalyxData3KO(11:100,6)).^2) + 10*sum((hz_10_3KO(101:111) - CalyxData3KO(101:111,2)).^2 + (hz_20_3KO(101:111) - CalyxData3KO(101:111,3)).^2 + (hz_50_3KO(101:111) - CalyxData3KO(101:111,4)).^2 + (hz_100_3KO(101:111) - CalyxData3KO(101:111,5)).^2 + (hz_200_3KO(101:111) - CalyxData3KO(101:111,6)).^2));


%DKO
C_7 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_DKO = state(end,:);

[~, ~, Fused_1_DKO] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_DKO = Fused_1_DKO/Fused_1_DKO(1);

[~, ~, Fused_10_DKO] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_DKO = Fused_10_DKO/Fused_10_DKO(1);

[~, ~, Fused_20_DKO] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_DKO = Fused_20_DKO/Fused_20_DKO(1);

[~, ~, Fused_50_DKO] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_DKO = Fused_50_DKO/Fused_50_DKO(1);

[~, ~, Fused_100_DKO] = stim_sim(stimulus_times_100, max_time_100, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_DKO = Fused_100_DKO/Fused_100_DKO(1);

[~, ~, Fused_200_DKO] = stim_sim(stimulus_times_200, max_time_200, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_DKO = Fused_200_DKO/Fused_200_DKO(1);

err_DKO = sqrt(2*sum((hz_1_DKO(1:10) - CalyxDataDKO(1:10,1)).^2 + (hz_10_DKO(1:10) - CalyxDataDKO(1:10,2)).^2 + (hz_20_DKO(1:10) - CalyxDataDKO(1:10,3)).^2 + (hz_50_DKO(1:10) - CalyxDataDKO(1:10,4)).^2 + (hz_100_DKO(1:10) - CalyxDataDKO(1:10,5)).^2 + (hz_200_DKO(1:10) - CalyxDataDKO(1:10,6)).^2) + sum((hz_1_DKO(11:100) - CalyxDataDKO(11:100,1)).^2 + (hz_10_DKO(11:100) - CalyxDataDKO(11:100,2)).^2 + (hz_20_DKO(11:100) - CalyxDataDKO(11:100,3)).^2 + (hz_50_DKO(11:100) - CalyxDataDKO(11:100,4)).^2 + (hz_100_DKO(11:100) - CalyxDataDKO(11:100,5)).^2 + (hz_200_DKO(11:100) - CalyxDataDKO(11:100,6)).^2) + 10*sum((hz_10_DKO(101:111) - CalyxDataDKO(101:111,2)).^2 + (hz_20_DKO(101:111) - CalyxDataDKO(101:111,3)).^2 + (hz_50_DKO(101:111) - CalyxDataDKO(101:111,4)).^2 + (hz_100_DKO(101:111) - CalyxDataDKO(101:111,5)).^2 + (hz_200_DKO(101:111) - CalyxDataDKO(101:111,6)).^2));

%7KO
C_3 = x(6);
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_7KO = state(end,:);

[~, ~, Fused_1_7KO] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_7KO, reserve_size, k_refill, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_7KO = Fused_1_7KO/Fused_1_7KO(1);

[~, ~, Fused_10_7KO] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_7KO, reserve_size, k_refill, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_7KO = Fused_10_7KO/Fused_10_7KO(1);

[~, ~, Fused_20_7KO] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_7KO, reserve_size, k_refill, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_7KO = Fused_20_7KO/Fused_20_7KO(1);

[~, ~, Fused_50_7KO] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_7KO, reserve_size, k_refill, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_7KO = Fused_50_7KO/Fused_50_7KO(1);

[~, ~, Fused_100_7KO] = stim_sim(stimulus_times_100, max_time_100, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_7KO, reserve_size, k_refill, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_7KO = Fused_100_7KO/Fused_100_7KO(1);

[~, ~, Fused_200_7KO] = stim_sim(stimulus_times_200, max_time_200, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_7KO, reserve_size, k_refill, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_7KO = Fused_200_7KO/Fused_200_7KO(1);

err_7KO = sqrt(2*sum((hz_1_7KO(1:10) - CalyxData7KO(1:10,1)).^2 + (hz_10_7KO(1:10) - CalyxData7KO(1:10,2)).^2 + (hz_20_7KO(1:10) - CalyxData7KO(1:10,3)).^2 + (hz_50_7KO(1:10) - CalyxData7KO(1:10,4)).^2 + (hz_100_7KO(1:10) - CalyxData7KO(1:10,5)).^2 + (hz_200_7KO(1:10) - CalyxData7KO(1:10,6)).^2) + sum((hz_1_7KO(11:100) - CalyxData7KO(11:100,1)).^2 + (hz_10_7KO(11:100) - CalyxData7KO(11:100,2)).^2 + (hz_20_7KO(11:100) - CalyxData7KO(11:100,3)).^2 + (hz_50_7KO(11:100) - CalyxData7KO(11:100,4)).^2 + (hz_100_7KO(11:100) - CalyxData7KO(11:100,5)).^2 + (hz_200_7KO(11:100) - CalyxData7KO(11:100,6)).^2) + 10*sum((hz_10_7KO(101:111) - CalyxData7KO(101:111,2)).^2 + (hz_20_7KO(101:111) - CalyxData7KO(101:111,3)).^2 + (hz_50_7KO(101:111) - CalyxData7KO(101:111,4)).^2 + (hz_100_7KO(101:111) - CalyxData7KO(101:111,5)).^2 + (hz_200_7KO(101:111) - CalyxData7KO(101:111,6)).^2));

err = (err_WT + err_3KO + err_7KO + err_DKO)/4;

cost = err + abs(err - err_WT)/25+ abs(err - err_3KO)/25+ abs(err - err_7KO)/25+ abs(err - err_DKO)/25;

disp(['Cost = ', num2str(cost), ', average error = ', num2str(err), ', WT error = ', num2str(err_WT), ', 3KO error = ', num2str(err_3KO), ', 7KO error = ', num2str(err_7KO), ', DKO error = ', num2str(err_DKO)])


function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS, reserve_size, k_refill, Syt3, Syt7)

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
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, Syt3(1,:), Syt7(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if ~isequal(stimulus_times, linspace(0,1000*99,100)) %1hz
       
        Recovery = [10 20 50 100 200 500 1000 2000 5000 10000 20000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [~,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, Syt3(i,:), Syt7(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,k_tether,k_untether,reserve_size,k_refill,Syt3,Syt7)

dydt(1,1) = -state(1)*(state(4)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1) - state(2)*k_tether*(1+Syt3+Syt7) + state(3)*k_untether;
dydt(3,1) = state(2)*k_tether*(1+Syt3+Syt7) - state(3)*k_untether;
dydt(4,1) = dydt(1,1) + (reserve_size-state(4))*k_refill;

end

function dydt = dState(t,state,k_docking,k_undocking,k_tether,k_untether,reserve_size,k_refill,Syt3,Syt7)

dydt(1,1) = -state(1)*(state(4)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1) - state(2)*k_tether*(1+Syt3(round(t/.1)+1)+Syt7(round(t/.1)+1)) + state(3)*k_untether;
dydt(3,1) = state(2)*k_tether*(1+Syt3(round(t/.1)+1)+Syt7(round(t/.1)+1)) - state(3)*k_untether;
dydt(4,1) = dydt(1,1) + (reserve_size-state(4))*k_refill;

end
end