CalyxDataWT = load('..\CalyxDataWT.mat').CalyxDataWT;
CalyxData3KO = load('..\CalyxData3KO.mat').CalyxData3KO;
CalyxData7KO = load('..\CalyxData7KO.mat').CalyxData7KO;
CalyxDataDKO = load('..\CalyxDataDKO.mat').CalyxDataDKO;
Calyx1HzCa = load('..\Calyx1HzCa.mat').Ca_sim;
Calyx10HzCa = load('..\Calyx10HzCa.mat').Ca_sim;
Calyx20HzCa = load('..\Calyx20HzCa.mat').Ca_sim;
Calyx50HzCa = load('..\Calyx50HzCa.mat').Ca_sim;
Calyx100HzCa = load('..\Calyx100HzCa.mat').Ca_sim;
Calyx200HzCa = load('..\Calyx200HzCa.mat').Ca_sim;
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

x = [0.66322,0.003346,0.00072235,8.687,5.9651e-05,0.063422,49.9845,10]; %best fit 2/15/21 cost=1.1365

p_release = x(1); 
k_docking = x(2);
k_undocking = x(3);
reserve_size = x(4);
k_refill = x(5);

Syt_pool_size = x(6);

C_3 = x(7);
C_7 = x(8);

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
state_0 = [1-Syt_pool_size; 0; Syt_pool_size; 0; reserve_size]; %[empty pool; bound pool; empty syt pool; bound syt pool; reserve pool]

%WT
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_WT = state(end,:);

[ts_1_WT, state_1_WT, Fused_1_WT] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_WT = Fused_1_WT/Fused_1_WT(1);

[ts_10_WT, state_10_WT, Fused_10_WT] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_WT = Fused_10_WT/Fused_10_WT(1);

[ts_20_WT, state_20_WT, Fused_20_WT] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_WT = Fused_20_WT/Fused_20_WT(1);

[ts_50_WT, state_50_WT, Fused_50_WT] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_WT = Fused_50_WT/Fused_50_WT(1);

[ts_100_WT, state_100_WT, Fused_100_WT] = stim_sim(stimulus_times_100, max_time_100, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_WT = Fused_100_WT/Fused_100_WT(1);

[ts_200_WT, state_200_WT, Fused_200_WT] = stim_sim(stimulus_times_200, max_time_200, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, Syt_pool_size, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_WT = Fused_200_WT/Fused_200_WT(1);

err_WT = sqrt(2*sum((hz_1_WT(1:10) - CalyxDataWT(1:10,1)).^2 + (hz_10_WT(1:10) - CalyxDataWT(1:10,2)).^2 + (hz_20_WT(1:10) - CalyxDataWT(1:10,3)).^2 + (hz_50_WT(1:10) - CalyxDataWT(1:10,4)).^2 + (hz_100_WT(1:10) - CalyxDataWT(1:10,5)).^2 + (hz_200_WT(1:10) - CalyxDataWT(1:10,6)).^2) + sum((hz_1_WT(11:100) - CalyxDataWT(11:100,1)).^2 + (hz_10_WT(11:100) - CalyxDataWT(11:100,2)).^2 + (hz_20_WT(11:100) - CalyxDataWT(11:100,3)).^2 + (hz_50_WT(11:100) - CalyxDataWT(11:100,4)).^2 + (hz_100_WT(11:100) - CalyxDataWT(11:100,5)).^2 + (hz_200_WT(11:100) - CalyxDataWT(11:100,6)).^2) + 10*sum((hz_10_WT(101:111) - CalyxDataWT(101:111,2)).^2 + (hz_20_WT(101:111) - CalyxDataWT(101:111,3)).^2 + (hz_50_WT(101:111) - CalyxDataWT(101:111,4)).^2 + (hz_100_WT(101:111) - CalyxDataWT(101:111,5)).^2 + (hz_200_WT(101:111) - CalyxDataWT(101:111,6)).^2));

%3KO
C_3 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_3KO = state(end,:);

[ts_1_3KO, state_1_3KO, Fused_1_3KO] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_3KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_3KO = Fused_1_3KO/Fused_1_3KO(1);

[ts_10_3KO, state_10_3KO, Fused_10_3KO] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_3KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_3KO = Fused_10_3KO/Fused_10_3KO(1);

[ts_20_3KO, state_20_3KO, Fused_20_3KO] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_3KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_3KO = Fused_20_3KO/Fused_20_3KO(1);

[ts_50_3KO, state_50_3KO, Fused_50_3KO] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_3KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_3KO = Fused_50_3KO/Fused_50_3KO(1);

[ts_100_3KO, state_100_3KO, Fused_100_3KO] = stim_sim(stimulus_times_100, max_time_100, p_release, k_docking, k_undocking, SS_3KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_3KO = Fused_100_3KO/Fused_100_3KO(1);

[ts_200_3KO, state_200_3KO, Fused_200_3KO] = stim_sim(stimulus_times_200, max_time_200, p_release, k_docking, k_undocking, SS_3KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_3KO = Fused_200_3KO/Fused_200_3KO(1);

err_3KO = sqrt(2*sum((hz_1_3KO(1:10) - CalyxData3KO(1:10,1)).^2 + (hz_10_3KO(1:10) - CalyxData3KO(1:10,2)).^2 + (hz_20_3KO(1:10) - CalyxData3KO(1:10,3)).^2 + (hz_50_3KO(1:10) - CalyxData3KO(1:10,4)).^2 + (hz_100_3KO(1:10) - CalyxData3KO(1:10,5)).^2 + (hz_200_3KO(1:10) - CalyxData3KO(1:10,6)).^2) + sum((hz_1_3KO(11:100) - CalyxData3KO(11:100,1)).^2 + (hz_10_3KO(11:100) - CalyxData3KO(11:100,2)).^2 + (hz_20_3KO(11:100) - CalyxData3KO(11:100,3)).^2 + (hz_50_3KO(11:100) - CalyxData3KO(11:100,4)).^2 + (hz_100_3KO(11:100) - CalyxData3KO(11:100,5)).^2 + (hz_200_3KO(11:100) - CalyxData3KO(11:100,6)).^2) + 10*sum((hz_10_3KO(101:111) - CalyxData3KO(101:111,2)).^2 + (hz_20_3KO(101:111) - CalyxData3KO(101:111,3)).^2 + (hz_50_3KO(101:111) - CalyxData3KO(101:111,4)).^2 + (hz_100_3KO(101:111) - CalyxData3KO(101:111,5)).^2 + (hz_200_3KO(101:111) - CalyxData3KO(101:111,6)).^2));

%DKO
C_7 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_DKO = state(end,:);

[ts_1_DKO, state_1_DKO, Fused_1_DKO] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_DKO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_DKO = Fused_1_DKO/Fused_1_DKO(1);

[ts_10_DKO, state_10_DKO, Fused_10_DKO] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_DKO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_DKO = Fused_10_DKO/Fused_10_DKO(1);

[ts_20_DKO, state_20_DKO, Fused_20_DKO] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_DKO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_DKO = Fused_20_DKO/Fused_20_DKO(1);

[ts_50_DKO, state_50_DKO, Fused_50_DKO] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_DKO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_DKO = Fused_50_DKO/Fused_50_DKO(1);

[ts_100_DKO, state_100_DKO, Fused_100_DKO] = stim_sim(stimulus_times_100, max_time_100, p_release, k_docking, k_undocking, SS_DKO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_DKO = Fused_100_DKO/Fused_100_DKO(1);

[ts_200_DKO, state_200_DKO, Fused_200_DKO] = stim_sim(stimulus_times_200, max_time_200, p_release, k_docking, k_undocking, SS_DKO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_DKO = Fused_200_DKO/Fused_200_DKO(1);

%7KO
C_3 = x(7);
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3(1), C_7*Calyx1HzSyt7(1)), [0 t_SS], state_0);

SS_7KO = state(end,:);

[ts_1_7KO, state_1_7KO, Fused_1_7KO] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_7KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx1HzSyt3, C_7*Calyx1HzSyt7);
hz_1_7KO = Fused_1_7KO/Fused_1_7KO(1);

[ts_10_7KO, state_10_7KO, Fused_10_7KO] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_7KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx10HzSyt3, C_7*Calyx10HzSyt7);
hz_10_7KO = Fused_10_7KO/Fused_10_7KO(1);

[ts_20_7KO, state_20_7KO, Fused_20_7KO] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_7KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx20HzSyt3, C_7*Calyx20HzSyt7);
hz_20_7KO = Fused_20_7KO/Fused_20_7KO(1);

[ts_50_7KO, state_50_7KO, Fused_50_7KO] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_7KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx50HzSyt3, C_7*Calyx50HzSyt7);
hz_50_7KO = Fused_50_7KO/Fused_50_7KO(1);

[ts_100_7KO, state_100_7KO, Fused_100_7KO] = stim_sim(stimulus_times_100, max_time_100, p_release, k_docking, k_undocking, SS_7KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx100HzSyt3, C_7*Calyx100HzSyt7);
hz_100_7KO = Fused_100_7KO/Fused_100_7KO(1);

[ts_200_7KO, state_200_7KO, Fused_200_7KO] = stim_sim(stimulus_times_200, max_time_200, p_release, k_docking, k_undocking, SS_7KO, reserve_size, k_refill, Syt_pool_size, C_3*Calyx200HzSyt3, C_7*Calyx200HzSyt7);
hz_200_7KO = Fused_200_7KO/Fused_200_7KO(1);

err_7KO = sqrt(2*sum((hz_1_7KO(1:10) - CalyxData7KO(1:10,1)).^2 + (hz_10_7KO(1:10) - CalyxData7KO(1:10,2)).^2 + (hz_20_7KO(1:10) - CalyxData7KO(1:10,3)).^2 + (hz_50_7KO(1:10) - CalyxData7KO(1:10,4)).^2 + (hz_100_7KO(1:10) - CalyxData7KO(1:10,5)).^2 + (hz_200_7KO(1:10) - CalyxData7KO(1:10,6)).^2) + sum((hz_1_7KO(11:100) - CalyxData7KO(11:100,1)).^2 + (hz_10_7KO(11:100) - CalyxData7KO(11:100,2)).^2 + (hz_20_7KO(11:100) - CalyxData7KO(11:100,3)).^2 + (hz_50_7KO(11:100) - CalyxData7KO(11:100,4)).^2 + (hz_100_7KO(11:100) - CalyxData7KO(11:100,5)).^2 + (hz_200_7KO(11:100) - CalyxData7KO(11:100,6)).^2) + 10*sum((hz_10_7KO(101:111) - CalyxData7KO(101:111,2)).^2 + (hz_20_7KO(101:111) - CalyxData7KO(101:111,3)).^2 + (hz_50_7KO(101:111) - CalyxData7KO(101:111,4)).^2 + (hz_100_7KO(101:111) - CalyxData7KO(101:111,5)).^2 + (hz_200_7KO(101:111) - CalyxData7KO(101:111,6)).^2));

err_DKO = sqrt(2*sum((hz_1_DKO(1:10) - CalyxDataDKO(1:10,1)).^2 + (hz_10_DKO(1:10) - CalyxDataDKO(1:10,2)).^2 + (hz_20_DKO(1:10) - CalyxDataDKO(1:10,3)).^2 + (hz_50_DKO(1:10) - CalyxDataDKO(1:10,4)).^2 + (hz_100_DKO(1:10) - CalyxDataDKO(1:10,5)).^2 + (hz_200_DKO(1:10) - CalyxDataDKO(1:10,6)).^2) + sum((hz_1_DKO(11:100) - CalyxDataDKO(11:100,1)).^2 + (hz_10_DKO(11:100) - CalyxDataDKO(11:100,2)).^2 + (hz_20_DKO(11:100) - CalyxDataDKO(11:100,3)).^2 + (hz_50_DKO(11:100) - CalyxDataDKO(11:100,4)).^2 + (hz_100_DKO(11:100) - CalyxDataDKO(11:100,5)).^2 + (hz_200_DKO(11:100) - CalyxDataDKO(11:100,6)).^2) + 10*sum((hz_10_DKO(101:111) - CalyxDataDKO(101:111,2)).^2 + (hz_20_DKO(101:111) - CalyxDataDKO(101:111,3)).^2 + (hz_50_DKO(101:111) - CalyxDataDKO(101:111,4)).^2 + (hz_100_DKO(101:111) - CalyxDataDKO(101:111,5)).^2 + (hz_200_DKO(101:111) - CalyxDataDKO(101:111,6)).^2));

err = (err_WT + err_3KO + err_7KO + err_DKO)/4;

cost = err + abs(err - err_WT)/25+ abs(err - err_3KO)/25+ abs(err - err_7KO)/25+ abs(err - err_DKO)/25;

disp(['Cost = ', num2str(cost), ', average error = ', num2str(err), ', WT error = ', num2str(err_WT), ', 3KO error = ', num2str(err_3KO), ', 7KO error = ', num2str(err_7KO), ', DKO error = ', num2str(err_DKO)])

%plot Ca and Syts
figure('Name','Ca & Syt3 Simulation (WT)','NumberTitle','off')
subplot(6,2,1)
semilogy(ts_1_WT,Calyx1HzCa(round(ts_1_WT/.1)+1))
title('1 Hz Ca')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

ax = gca;
ax.XRuler.Exponent = 0;

subplot(6,2,2)
plot(ts_1_WT,Calyx1HzSyt3(round(ts_1_WT/.1)+1))
hold on
plot(ts_1_WT,Calyx1HzSyt7(round(ts_1_WT/.1)+1))
title('1 Hz Syt')
xlabel('time (ms)')
ylabel('Bound Syt')

ax = gca;
ax.XRuler.Exponent = 0;

subplot(6,2,3)
semilogy(ts_10_WT,Calyx10HzCa(1,round(ts_10_WT/.1)+1))
title('10 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(6,2,4)
plot(ts_10_WT,Calyx10HzSyt3(1,round(ts_10_WT/.1)+1))
hold on
plot(ts_10_WT,Calyx10HzSyt7(1,round(ts_10_WT/.1)+1))
title('10 Hz')
xlabel('time (ms)')
ylabel('Bound Syt')

subplot(6,2,5)
semilogy(ts_20_WT,Calyx20HzCa(1,round(ts_20_WT/.1)+1))
title('20 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(6,2,6)
plot(ts_20_WT,Calyx20HzSyt3(1,round(ts_20_WT/.1)+1))
hold on
plot(ts_20_WT,Calyx20HzSyt7(1,round(ts_20_WT/.1)+1))
title('20 Hz')
xlabel('time (ms)')
ylabel('Bound Syt')

subplot(6,2,7)
semilogy(ts_50_WT,Calyx50HzCa(1,round(ts_50_WT/.1)+1))
title('50 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(6,2,8)
plot(ts_50_WT,Calyx50HzSyt3(1,round(ts_50_WT/.1)+1))
hold on
plot(ts_50_WT,Calyx50HzSyt7(1,round(ts_50_WT/.1)+1))
title('50 Hz')
xlabel('time (ms)')
ylabel('Bound Syt')

subplot(6,2,9)
semilogy(ts_100_WT,Calyx100HzCa(1,round(ts_100_WT/.1)+1))
title('100 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(6,2,10)
plot(ts_100_WT,Calyx100HzSyt3(1,round(ts_100_WT/.1)+1))
hold on
plot(ts_100_WT,Calyx100HzSyt7(1,round(ts_100_WT/.1)+1))
title('100 Hz')
xlabel('time (ms)')
ylabel('Bound Syt')

subplot(6,2,11)
semilogy(ts_200_WT,Calyx200HzCa(1,round(ts_200_WT/.1)+1))
title('200 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(6,2,12)
plot(ts_200_WT,Calyx200HzSyt3(1,round(ts_200_WT/.1)+1))
hold on
plot(ts_200_WT,Calyx200HzSyt7(1,round(ts_200_WT/.1)+1))
title('200 Hz')
xlabel('time (ms)')
ylabel('Bound Syt')

%state plots
figure('Name','State Simulations','NumberTitle','off')
subplot(6,4,1)
plot(ts_1_WT,state_1_WT(:,1),'color',[255, 165, 0]/255)
title('WT 1 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_1_WT,state_1_WT(:,2),'-b')
plot(ts_1_WT,state_1_WT(:,3),'-k')

ax = gca;
ax.XRuler.Exponent = 0;

subplot(6,4,5)
plot(ts_10_WT,state_10_WT(:,1),'color',[255, 165, 0]/255)
title('WT 10 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_10_WT,state_10_WT(:,2),'-b')
plot(ts_10_WT,state_10_WT(:,3),'-k')

subplot(6,4,9)
plot(ts_20_WT,state_20_WT(:,1),'color',[255, 165, 0]/255)
title('WT 20 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_20_WT,state_20_WT(:,2),'-b')
plot(ts_20_WT,state_20_WT(:,3),'-k')

subplot(6,4,13)
plot(ts_50_WT,state_50_WT(:,1),'color',[255, 165, 0]/255)
title('WT 50 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_50_WT,state_50_WT(:,2),'-b')
plot(ts_50_WT,state_50_WT(:,3),'-k')

subplot(6,4,17)
plot(ts_100_WT,state_100_WT(:,1),'color',[255, 165, 0]/255)
title('WT 100 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_100_WT,state_100_WT(:,2),'-b')
plot(ts_100_WT,state_100_WT(:,3),'-k')

subplot(6,4,21)
plot(ts_200_WT,state_200_WT(:,1),'color',[255, 165, 0]/255)
title('WT 200 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_200_WT,state_200_WT(:,2),'-b')
plot(ts_200_WT,state_200_WT(:,3),'-k')

subplot(6,4,2)
plot(ts_1_3KO,state_1_3KO(:,1),'color',[255, 165, 0]/255)
title('3KO 1 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_1_3KO,state_1_3KO(:,2),'-b')
plot(ts_1_3KO,state_1_3KO(:,3),'-k')

ax = gca;
ax.XRuler.Exponent = 0;

legend({'Empty Sites','Filled Sites','Reserve Vesicles'},'Location','Best')

subplot(6,4,6)
plot(ts_10_3KO,state_10_3KO(:,1),'color',[255, 165, 0]/255)
title('3KO 10 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_10_3KO,state_10_3KO(:,2),'-b')
plot(ts_10_3KO,state_10_3KO(:,3),'-k')

subplot(6,4,10)
plot(ts_20_3KO,state_20_3KO(:,1),'color',[255, 165, 0]/255)
title('3KO 20 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_20_3KO,state_20_3KO(:,2),'-b')
plot(ts_20_3KO,state_20_3KO(:,3),'-k')

subplot(6,4,14)
plot(ts_50_3KO,state_50_3KO(:,1),'color',[255, 165, 0]/255)
title('3KO 50 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_50_3KO,state_50_3KO(:,2),'-b')
plot(ts_50_3KO,state_50_3KO(:,3),'-k')

subplot(6,4,18)
plot(ts_100_3KO,state_100_3KO(:,1),'color',[255, 165, 0]/255)
title('3KO 100 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_100_3KO,state_100_3KO(:,2),'-b')
plot(ts_100_3KO,state_100_3KO(:,3),'-k')

subplot(6,4,22)
plot(ts_200_3KO,state_200_3KO(:,1),'color',[255, 165, 0]/255)
title('3KO 200 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_200_3KO,state_200_3KO(:,2),'-b')
plot(ts_200_3KO,state_200_3KO(:,3),'-k')

subplot(6,4,3)
plot(ts_1_7KO,state_1_7KO(:,1),'color',[255, 165, 0]/255)
title('7KO 1 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_1_7KO,state_1_7KO(:,2),'-b')
plot(ts_1_7KO,state_1_7KO(:,3),'-k')

ax = gca;
ax.XRuler.Exponent = 0;

subplot(6,4,7)
plot(ts_10_7KO,state_10_7KO(:,1),'color',[255, 165, 0]/255)
title('7KO 10 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_10_7KO,state_10_7KO(:,2),'-b')
plot(ts_10_7KO,state_10_7KO(:,3),'-k')

subplot(6,4,11)
plot(ts_20_7KO,state_20_7KO(:,1),'color',[255, 165, 0]/255)
title('7KO 20 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_20_7KO,state_20_7KO(:,2),'-b')
plot(ts_20_7KO,state_20_7KO(:,3),'-k')

subplot(6,4,15)
plot(ts_50_7KO,state_50_7KO(:,1),'color',[255, 165, 0]/255)
title('7KO 50 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_50_7KO,state_50_7KO(:,2),'-b')
plot(ts_50_7KO,state_50_7KO(:,3),'-k')

subplot(6,4,19)
plot(ts_100_7KO,state_100_7KO(:,1),'color',[255, 165, 0]/255)
title('7KO 100 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_100_7KO,state_100_7KO(:,2),'-b')
plot(ts_100_7KO,state_100_7KO(:,3),'-k')

subplot(6,4,23)
plot(ts_200_7KO,state_200_7KO(:,1),'color',[255, 165, 0]/255)
title('7KO 200 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_200_7KO,state_200_7KO(:,2),'-b')
plot(ts_200_7KO,state_200_7KO(:,3),'-k')

subplot(6,4,4)
plot(ts_1_DKO,state_1_DKO(:,1),'color',[255, 165, 0]/255)
title('DKO 1 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_1_DKO,state_1_DKO(:,2),'-b')
plot(ts_1_DKO,state_1_DKO(:,3),'-k')

ax = gca;
ax.XRuler.Exponent = 0;

subplot(6,4,8)
plot(ts_10_DKO,state_10_DKO(:,1),'color',[255, 165, 0]/255)
title('DKO 10 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_10_DKO,state_10_DKO(:,2),'-b')
plot(ts_10_DKO,state_10_DKO(:,3),'-k')

subplot(6,4,12)
plot(ts_20_DKO,state_20_DKO(:,1),'color',[255, 165, 0]/255)
title('DKO 20 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_20_DKO,state_20_DKO(:,2),'-b')
plot(ts_20_DKO,state_20_DKO(:,3),'-k')

subplot(6,4,16)
plot(ts_50_DKO,state_50_DKO(:,1),'color',[255, 165, 0]/255)
title('DKO 50 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_50_DKO,state_50_DKO(:,2),'-b')
plot(ts_50_DKO,state_50_DKO(:,3),'-k')

subplot(6,4,20)
plot(ts_100_DKO,state_100_DKO(:,1),'color',[255, 165, 0]/255)
title('DKO 100 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_100_DKO,state_100_DKO(:,2),'-b')
plot(ts_100_DKO,state_100_DKO(:,3),'-k')

subplot(6,4,24)
plot(ts_200_DKO,state_200_DKO(:,1),'color',[255, 165, 0]/255)
title('DKO 200 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_200_DKO,state_200_DKO(:,2),'-b')
plot(ts_200_DKO,state_200_DKO(:,3),'-k')


%plot simulated data
rec = [10 20 50 100 200 500 1000 2000 5000 10000 20000];

figure('Name','WT Simulated vs Collected Data','NumberTitle','off')
subplot(1,2,1)
plot(CalyxDataWT(1:100,1),'-k')
title('WT Trains')
xlabel('Stim #')
ylabel('EPSC (norm)')
set(gca,'xlim',[1 100])
set(gca,'ylim',[0 1.2])
hold on
plot(CalyxDataWT(1:100,2),'-k')
plot(CalyxDataWT(1:100,3),'-k')
plot(CalyxDataWT(1:100,4),'-k')
plot(CalyxDataWT(1:100,5),'-k')
plot(CalyxDataWT(1:100,6),'-k')
plot(hz_1_WT,'-','color',[128, 128, 128]/255)
plot(hz_10_WT(1:100),'-','color',[128, 128, 128]/255)
plot(hz_20_WT(1:100),'-','color',[128, 128, 128]/255)
plot(hz_50_WT(1:100),'-','color',[128, 128, 128]/255)
plot(hz_100_WT(1:100),'-','color',[128, 128, 128]/255)
plot(hz_200_WT(1:100),'-','color',[128, 128, 128]/255)

subplot(1,2,2)
title('WT Recovery')
xlabel('t(ms)')
ylabel('EPSC (norm)')
semilogx(rec,CalyxDataWT(101:111,2),'-','color',[0, 0, 0]/255)
hold on
set(gca,'xlim',[10 25000])
set(gca,'ylim',[0 1.2])
semilogx(rec,CalyxDataWT(101:111,3),'-','color',[255, 0, 0]/255)
semilogx(rec,CalyxDataWT(101:111,4),'-','color',[0, 0, 255]/255)
semilogx(rec,CalyxDataWT(101:111,5),'-','color',[255, 128, 0]/255)
semilogx(rec,CalyxDataWT(101:111,6),'-','color',[255, 0, 255]/255)
semilogx(rec,hz_10_WT(101:111),'o','color',[255, 0, 255]/255)
semilogx(rec,hz_20_WT(101:111),'o','color',[255, 0, 0]/255)
semilogx(rec,hz_50_WT(101:111),'o','color',[0, 0, 0]/255)
semilogx(rec,hz_100_WT(101:111),'o','color',[255, 128, 0]/255)
semilogx(rec,hz_200_WT(101:111),'o','color',[255, 0, 255]/255)
legend({'10 Hz Data', '20 Hz Data', '50 Hz Data', '100 Hz Data', '200 Hz Data'},'Location','Best')

figure('Name','3KO Simulated vs Collected Data','NumberTitle','off')
subplot(1,2,1)
plot(CalyxData3KO(1:100,1),'-k')
title('3KO Trains')
xlabel('Stim #')
ylabel('EPSC (norm)')
set(gca,'xlim',[1 100])
set(gca,'ylim',[0 1.2])
hold on
plot(CalyxData3KO(1:100,2),'-k')
plot(CalyxData3KO(1:100,3),'-k')
plot(CalyxData3KO(1:100,4),'-k')
plot(CalyxData3KO(1:100,5),'-k')
plot(CalyxData3KO(1:100,6),'-k')
plot(hz_1_3KO,'-','color',[128, 128, 128]/255)
plot(hz_10_3KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_20_3KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_50_3KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_100_3KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_200_3KO(1:100),'-','color',[128, 128, 128]/255)

subplot(1,2,2)
title('3KO Recovery')
xlabel('t(ms)')
ylabel('EPSC (norm)')
semilogx(rec,CalyxData3KO(101:111,2),'-','color',[0, 0, 0]/255)
hold on
set(gca,'xlim',[10 25000])
set(gca,'ylim',[0 1.2])
semilogx(rec,CalyxData3KO(101:111,3),'-','color',[255, 0, 0]/255)
semilogx(rec,CalyxData3KO(101:111,4),'-','color',[0, 0, 255]/255)
semilogx(rec,CalyxData3KO(101:111,5),'-','color',[255, 128, 0]/255)
semilogx(rec,CalyxData3KO(101:111,6),'-','color',[255, 0, 255]/255)
semilogx(rec,hz_10_3KO(101:111),'o','color',[255, 0, 255]/255)
semilogx(rec,hz_20_3KO(101:111),'o','color',[255, 0, 0]/255)
semilogx(rec,hz_50_3KO(101:111),'o','color',[0, 0, 0]/255)
semilogx(rec,hz_100_3KO(101:111),'o','color',[255, 128, 0]/255)
semilogx(rec,hz_200_3KO(101:111),'o','color',[255, 0, 255]/255)
legend({'10 Hz Data', '20 Hz Data', '50 Hz Data', '100 Hz Data', '200 Hz Data'},'Location','Best')

figure('Name','7KO Simulated vs Collected Data','NumberTitle','off')
subplot(1,2,1)
plot(CalyxData7KO(1:100,1),'-k')
title('7KO Trains')
xlabel('Stim #')
ylabel('EPSC (norm)')
set(gca,'xlim',[1 100])
set(gca,'ylim',[0 1.2])
hold on
plot(CalyxData7KO(1:100,2),'-k')
plot(CalyxData7KO(1:100,3),'-k')
plot(CalyxData7KO(1:100,4),'-k')
plot(CalyxData7KO(1:100,5),'-k')
plot(CalyxData7KO(1:100,6),'-k')
plot(hz_1_7KO,'-','color',[128, 128, 128]/255)
plot(hz_10_7KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_20_7KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_50_7KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_100_7KO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_200_7KO(1:100),'-','color',[128, 128, 128]/255)

subplot(1,2,2)
title('7KO Recovery')
xlabel('t(ms)')
ylabel('EPSC (norm)')
semilogx(rec,CalyxData7KO(101:111,2),'-','color',[0, 0, 0]/255)
hold on
set(gca,'xlim',[10 25000])
set(gca,'ylim',[0 1.2])
semilogx(rec,CalyxData7KO(101:111,3),'-','color',[255, 0, 0]/255)
semilogx(rec,CalyxData7KO(101:111,4),'-','color',[0, 0, 255]/255)
semilogx(rec,CalyxData7KO(101:111,5),'-','color',[255, 128, 0]/255)
semilogx(rec,CalyxData7KO(101:111,6),'-','color',[255, 0, 255]/255)
semilogx(rec,hz_10_7KO(101:111),'o','color',[255, 0, 255]/255)
semilogx(rec,hz_20_7KO(101:111),'o','color',[255, 0, 0]/255)
semilogx(rec,hz_50_7KO(101:111),'o','color',[0, 0, 0]/255)
semilogx(rec,hz_100_7KO(101:111),'o','color',[255, 128, 0]/255)
semilogx(rec,hz_200_7KO(101:111),'o','color',[255, 0, 255]/255)
legend({'10 Hz Data', '20 Hz Data', '50 Hz Data', '100 Hz Data', '200 Hz Data'},'Location','Best')

figure('Name','DKO Simulated vs Collected Data','NumberTitle','off')
subplot(1,2,1)
plot(CalyxDataDKO(1:100,1),'-k')
title('DKO Trains')
xlabel('Stim #')
ylabel('EPSC (norm)')
set(gca,'xlim',[1 100])
set(gca,'ylim',[0 1.2])
hold on
plot(CalyxDataDKO(1:100,2),'-k')
plot(CalyxDataDKO(1:100,3),'-k')
plot(CalyxDataDKO(1:100,4),'-k')
plot(CalyxDataDKO(1:100,5),'-k')
plot(CalyxDataDKO(1:100,6),'-k')
plot(hz_1_DKO,'-','color',[128, 128, 128]/255)
plot(hz_10_DKO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_20_DKO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_50_DKO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_100_DKO(1:100),'-','color',[128, 128, 128]/255)
plot(hz_200_DKO(1:100),'-','color',[128, 128, 128]/255)

subplot(1,2,2)
title('DKO Recovery')
xlabel('t(ms)')
ylabel('EPSC (norm)')
semilogx(rec,CalyxDataDKO(101:111,2),'-','color',[0, 0, 0]/255)
hold on
set(gca,'xlim',[10 25000])
set(gca,'ylim',[0 1.2])
semilogx(rec,CalyxDataDKO(101:111,3),'-','color',[255, 0, 0]/255)
semilogx(rec,CalyxDataDKO(101:111,4),'-','color',[0, 0, 255]/255)
semilogx(rec,CalyxDataDKO(101:111,5),'-','color',[255, 128, 0]/255)
semilogx(rec,CalyxDataDKO(101:111,6),'-','color',[255, 0, 255]/255)
semilogx(rec,hz_10_DKO(101:111),'o','color',[255, 0, 255]/255)
semilogx(rec,hz_20_DKO(101:111),'o','color',[255, 0, 0]/255)
semilogx(rec,hz_50_DKO(101:111),'o','color',[0, 0, 0]/255)
semilogx(rec,hz_100_DKO(101:111),'o','color',[255, 128, 0]/255)
semilogx(rec,hz_200_DKO(101:111),'o','color',[255, 0, 255]/255)
legend({'10 Hz Data', '20 Hz Data', '50 Hz Data', '100 Hz Data', '200 Hz Data'},'Location','Best')


function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, SS, reserve_size, k_refill, Syt_pool_size, Syt3, Syt7)

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
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, Syt3(1,:), Syt7(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [~,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt_pool_size, Syt3(i,:), Syt7(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release + pre_stim(4)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,reserve_size,k_refill,Syt_pool_size,Syt3,Syt7)

dydt(1,1) = -(state(1)/(1-Syt_pool_size))*(state(5)/reserve_size)*k_docking + (state(2)/(1-Syt_pool_size))*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = -(state(3)/Syt_pool_size)*(state(5)/reserve_size)*k_docking*Syt3*Syt7 + (state(4)/Syt_pool_size)*k_undocking;
dydt(4,1) = -dydt(3,1);
dydt(5,1) = dydt(1,1) + dydt(3,1) + (reserve_size-state(5))*k_refill;

end

function dydt = dState(t,state,k_docking,k_undocking,reserve_size,k_refill,Syt_pool_size,Syt3,Syt7)

dydt(1,1) = -(state(1)/(1-Syt_pool_size))*(state(5)/reserve_size)*k_docking + (state(2)/(1-Syt_pool_size))*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = -(state(3)/Syt_pool_size)*(state(5)/reserve_size)*k_docking*Syt3(round(t/.1)+1)*Syt7(round(t/.1)+1) + (state(4)/Syt_pool_size)*k_undocking;
dydt(4,1) = -dydt(3,1);
dydt(5,1) = dydt(1,1) + dydt(3,1) + (reserve_size-state(5))*k_refill;

end