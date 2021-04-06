CalyxDataWT = load('..\CalyxDataWT.mat').CalyxDataWT;
CalyxData3KO = load('..\CalyxData3KO.mat').CalyxData3KO;
CalyxData7KO = load('..\CalyxData7KO.mat').CalyxData7KO;
CalyxDataDKO = load('..\CalyxDataDKO.mat').CalyxDataDKO;

%x = [0.2807,0.5991,0.0027,0.0156,0.7505,0.1754,43.1915,0.0338]; %cost = 8.9881
%x = [0.99931,0.92675,0.00038824,0.098077,0.0095092,0.0019097,44.6718,0.13843]; %with an error of 2.9304
%x = [0.42148,0.77812,0.0080469,0.0001,0.0002,5e-05,41.4844,1e-05]; %with an error of 2.5059
x = [0.95231,0.86974,0.0097268,0.098323,0.0049815,0.00068531,199.8649,0.30722]; %err = 2.1805

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
C_3 = 0;
C_7 = 0;
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill), [0 t_SS], state_0);

SS_DKO = state(end,:);

[ts_1_DKO, state_1_DKO, Fused_1_DKO] = stim_sim(stimulus_times_1, max_time_1, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_1_DKO = Fused_1_DKO/Fused_1_DKO(1);

[ts_10_DKO, state_10_DKO, Fused_10_DKO] = stim_sim(stimulus_times_10, max_time_10, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_10_DKO = Fused_10_DKO/Fused_10_DKO(1);

[ts_20_DKO, state_20_DKO, Fused_20_DKO] = stim_sim(stimulus_times_20, max_time_20, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_20_DKO = Fused_20_DKO/Fused_20_DKO(1);

[ts_50_DKO, state_50_DKO, Fused_50_DKO] = stim_sim(stimulus_times_50, max_time_50, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_50_DKO = Fused_50_DKO/Fused_50_DKO(1);

[ts_100_DKO, state_100_DKO, Fused_100_DKO] = stim_sim(stimulus_times_100, max_time_100, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_100_DKO = Fused_100_DKO/Fused_100_DKO(1);

[ts_200_DKO, state_200_DKO, Fused_200_DKO] = stim_sim(stimulus_times_200, max_time_200, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS_DKO, reserve_size, k_refill);
hz_200_DKO = Fused_200_DKO/Fused_200_DKO(1);

err_DKO = sqrt(2*sum((hz_1_DKO(1:10) - CalyxDataDKO(1:10,1)).^2 + (hz_10_DKO(1:10) - CalyxDataDKO(1:10,2)).^2 + (hz_20_DKO(1:10) - CalyxDataDKO(1:10,3)).^2 + (hz_50_DKO(1:10) - CalyxDataDKO(1:10,4)).^2 + (hz_100_DKO(1:10) - CalyxDataDKO(1:10,5)).^2 + (hz_200_DKO(1:10) - CalyxDataDKO(1:10,6)).^2) + sum((hz_1_DKO(11:100) - CalyxDataDKO(11:100,1)).^2 + (hz_10_DKO(11:100) - CalyxDataDKO(11:100,2)).^2 + (hz_20_DKO(11:100) - CalyxDataDKO(11:100,3)).^2 + (hz_50_DKO(11:100) - CalyxDataDKO(11:100,4)).^2 + (hz_100_DKO(11:100) - CalyxDataDKO(11:100,5)).^2 + (hz_200_DKO(11:100) - CalyxDataDKO(11:100,6)).^2) + 10*sum((hz_10_DKO(101:111) - CalyxDataDKO(101:111,2)).^2 + (hz_20_DKO(101:111) - CalyxDataDKO(101:111,3)).^2 + (hz_50_DKO(101:111) - CalyxDataDKO(101:111,4)).^2 + (hz_100_DKO(101:111) - CalyxDataDKO(101:111,5)).^2 + (hz_200_DKO(101:111) - CalyxDataDKO(101:111,6)).^2));

err = err_DKO;

cost = err + abs(err - err_DKO)/5;

disp(['Cost = ', num2str(cost), ', DKO error = ', num2str(err_DKO)])

%plot simulated data
rec = [10 20 50 100 200 500 1000 2000 5000 10000 20000];

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
% hold on
% set(gca,'xlim',[10 25000])
% set(gca,'ylim',[0 1.2])
% % semilogx(rec,hz_10_DKO(101:111),'o','color',[255, 0, 255]/255)
% % semilogx(rec,hz_20_DKO(101:111),'o','color',[255, 0, 0]/255)
% % semilogx(rec,hz_50_DKO(101:111),'o','color',[0, 0, 0]/255)
% % semilogx(rec,hz_100_DKO(101:111),'o','color',[255, 128, 0]/255)
semilogx(rec,hz_200_DKO(101:111),'o','color',[255, 0, 255]/255)
legend({'10 Hz Data', '20 Hz Data', '50 Hz Data', '100 Hz Data', '200 Hz Data'},'Location','Best')
% legend({'200 Hz Data'},'Location','Best')

function [ts, state, Fused] = stim_sim(stimulus_times, ~, p_release_dock, p_release_tether, k_docking, k_undocking, k_tether, k_untether, SS, reserve_size, k_refill)

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