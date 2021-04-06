CFDataWT = load('..\CFDataWT.mat').CFDataWT;
CFDataKO = load('..\CFData3KO.mat').CFData3KO;
CF1HzCa = load('..\CF1HzCa.mat').Ca_sim;
CF10HzCa = load('..\CF10HzCa.mat').Ca_sim;
CF20HzCa = load('..\CF20HzCa.mat').Ca_sim;
CF50HzCa = load('..\CF50HzCa.mat').Ca_sim;
CF1HzSyt3 = load('..\CF1HzSyt3.mat').Syt3;
CF10HzSyt3 = load('..\CF10HzSyt3.mat').Syt3;
CF20HzSyt3 = load('..\CF20HzSyt3.mat').Syt3;
CF50HzSyt3 = load('..\CF50HzSyt3.mat').Syt3;

plotCaSyts = 0;
plotStates = 0;
%x = [0.70303,0.0033015,0.00010035,6.7483,0.0001,0]; %KO best fit 1/13/2021 err = 0.69845
x = [0.737502,0.002,0.00146176,14.2557,3.5552e-05,6.48438]; %best fit 2/28/21 cost = 0.91255

p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);
k_refill = x(5);

C_3 = x(6); %multiplicative factor of syt3 effect

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

%WT
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3*CF1HzSyt3(1)), [0 t_SS], state_0);

SS_WT = state(end,:);

[ts_1_WT, state_1_WT, Fused_1_WT] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, C_3*CF1HzSyt3);
hz_1_WT = Fused_1_WT/Fused_1_WT(1);

[ts_10_WT, state_10_WT, Fused_10_WT] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, C_3*CF10HzSyt3);
hz_10_WT = Fused_10_WT/Fused_10_WT(1);

[ts_20_WT, state_20_WT, Fused_20_WT] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, C_3*CF20HzSyt3);
hz_20_WT = Fused_20_WT/Fused_20_WT(1);

[ts_50_WT, state_50_WT, Fused_50_WT] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill, C_3*CF50HzSyt3);
hz_50_WT = Fused_50_WT/Fused_50_WT(1);

err_WT = sqrt(2*sum((hz_1_WT(1:5) - CFDataWT(1:5,1)).^2 + (hz_10_WT(1:5) - CFDataWT(1:5,2)).^2 + (hz_20_WT(1:5) - CFDataWT(1:5,3)).^2 + (hz_50_WT(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_1_WT(6:20) - CFDataWT(6:20,1)).^2 + (hz_10_WT(6:20) - CFDataWT(6:20,2)).^2 + (hz_20_WT(6:20) - CFDataWT(6:20,3)).^2 + (hz_50_WT(6:20) - CFDataWT(6:20,4)).^2) + 10*sum((hz_50_WT(21:30) - CFDataWT(21:30,4)).^2));

%KO
C_3 = 0;
[t0,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3*CF1HzSyt3(1)), [0 t_SS], state_0);

SS_KO = state(end,:);

[ts_1_KO, state_1_KO, Fused_1_KO] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, C_3*CF1HzSyt3);
hz_1_KO = Fused_1_KO/Fused_1_KO(1);

[ts_10_KO, state_10_KO, Fused_10_KO] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, C_3*CF10HzSyt3);
hz_10_KO = Fused_10_KO/Fused_10_KO(1);

[ts_20_KO, state_20_KO, Fused_20_KO] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, C_3*CF20HzSyt3);
hz_20_KO = Fused_20_KO/Fused_20_KO(1);

[ts_50_KO, state_50_KO, Fused_50_KO] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_KO, reserve_size, k_refill, C_3*CF50HzSyt3);
hz_50_KO = Fused_50_KO/Fused_50_KO(1);

err_KO = sqrt(2*sum((hz_1_KO(1:5) - CFDataKO(1:5,1)).^2 + (hz_10_KO(1:5) - CFDataKO(1:5,2)).^2 + (hz_20_KO(1:5) - CFDataKO(1:5,3)).^2 + (hz_50_KO(1:5) - CFDataKO(1:5,4)).^2) + sum((hz_1_KO(6:20) - CFDataKO(6:20,1)).^2 + (hz_10_KO(6:20) - CFDataKO(6:20,2)).^2 + (hz_20_KO(6:20) - CFDataKO(6:20,3)).^2 + (hz_50_KO(6:20) - CFDataKO(6:20,4)).^2) + 10*sum((hz_50_KO(21:30) - CFDataKO(21:30,4)).^2));

err = (err_WT + err_KO)/2;

cost = err + abs(err_WT - err_KO)/10;

disp(['Cost = ', num2str(cost), ', average error = ', num2str(err), ', WT error = ', num2str(err_WT), ', KO error = ', num2str(err_KO)])

if plotCaSyts == 1
    %plot Ca and Syts
    figure('Name','Ca & Syt3 Simulation (WT)','NumberTitle','off')
    subplot(4,2,1)
    semilogy(ts_1_WT,CF1HzCa(round(ts_1_WT/.1)+1))
    title('1 Hz Ca')
    xlabel('time (ms)')
    ylabel('Ca Conc. (M)')

    ax = gca;
    ax.XRuler.Exponent = 0;

    subplot(4,2,2)
    plot(ts_1_WT,CF1HzSyt3(round(ts_1_WT/.1)+1))
    title('1 Hz Syt3')
    xlabel('time (ms)')
    ylabel('Bound Syt3')

    ax = gca;
    ax.XRuler.Exponent = 0;

    subplot(4,2,3)
    semilogy(ts_10_WT,CF10HzCa(round(ts_10_WT/.1)+1))
    title('10 Hz')
    xlabel('time (ms)')
    ylabel('Ca Conc. (M)')

    subplot(4,2,4)
    plot(ts_10_WT,CF10HzSyt3(round(ts_10_WT/.1)+1))
    title('10 Hz')
    xlabel('time (ms)')
    ylabel('Bound Syt3')

    subplot(4,2,5)
    semilogy(ts_20_WT,CF20HzCa(round(ts_20_WT/.1)+1))
    title('20 Hz')
    xlabel('time (ms)')
    ylabel('Ca Conc. (M)')

    subplot(4,2,6)
    plot(ts_20_WT,CF20HzSyt3(round(ts_20_WT/.1)+1))
    title('20 Hz')
    xlabel('time (ms)')
    ylabel('Bound Syt3')

    subplot(4,2,7)
    semilogy(ts_50_WT,CF50HzCa(1,round(ts_50_WT/.1)+1))
    title('50 Hz')
    xlabel('time (ms)')
    ylabel('Ca Conc. (M)')

    subplot(4,2,8)
    plot(ts_50_WT,CF50HzSyt3(1,round(ts_50_WT/.1)+1))
    title('50 Hz')
    xlabel('time (ms)')
    ylabel('Bound Syt3')
end

if plotStates  == 1
    %state plots
    figure('Name','State Simulations','NumberTitle','off')
    subplot(4,2,1)
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

    subplot(4,2,3)
    plot(ts_10_WT,state_10_WT(:,1),'color',[255, 165, 0]/255)
    title('WT 10 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_10_WT,state_10_WT(:,2),'-b')
    plot(ts_10_WT,state_10_WT(:,3),'-k')

    subplot(4,2,5)
    plot(ts_20_WT,state_20_WT(:,1),'color',[255, 165, 0]/255)
    title('WT 20 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_20_WT,state_20_WT(:,2),'-b')
    plot(ts_20_WT,state_20_WT(:,3),'-k')

    subplot(4,2,7)
    plot(ts_50_WT,state_50_WT(:,1),'color',[255, 165, 0]/255)
    title('WT 50 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_50_WT,state_50_WT(:,2),'-b')
    plot(ts_50_WT,state_50_WT(:,3),'-k')

    subplot(4,2,2)
    plot(ts_1_KO,state_1_KO(:,1),'color',[255, 165, 0]/255)
    title('KO 1 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_1_KO,state_1_KO(:,2),'-b')
    plot(ts_1_KO,state_1_KO(:,3),'-k')

    ax = gca;
    ax.XRuler.Exponent = 0;

    legend({'Empty Sites','Filled Sites','Reserve Vesicles'},'Location','Best')

    subplot(4,2,4)
    plot(ts_10_KO,state_10_KO(:,1),'color',[255, 165, 0]/255)
    title('KO 10 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_10_KO,state_10_KO(:,2),'-b')
    plot(ts_10_KO,state_10_KO(:,3),'-k')

    subplot(4,2,6)
    plot(ts_20_KO,state_20_KO(:,1),'color',[255, 165, 0]/255)
    title('KO 20 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_20_KO,state_20_KO(:,2),'-b')
    plot(ts_20_KO,state_20_KO(:,3),'-k')

    subplot(4,2,8)
    plot(ts_50_KO,state_50_KO(:,1),'color',[255, 165, 0]/255)
    title('KO 50 Hz')
    xlabel('time (ms)')
    ylabel('Value')
    set(gca,'ylim',[0 reserve_size])
    hold on
    plot(ts_50_KO,state_50_KO(:,2),'-b')
    plot(ts_50_KO,state_50_KO(:,3),'-k')
end

%plot simulated data
rec = [50 100 200 350 500 750 1000 2000 5000 10000];

figure('Name','Simulated vs Collected Data','NumberTitle','off')
subplot(5,2,1)
plot(CFDataWT(1:20,1),'-k')
title('WT 1 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_1_WT,'ko','Markersize',5)
legend({'Data','Simulation'},'Location','Best')

subplot(5,2,3)
plot(CFDataWT(1:20,2),'-k')
title('WT 10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10_WT,'ko','Markersize',5)

subplot(5,2,5)
plot(CFDataWT(1:20,3),'-k')
title('WT 20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20_WT,'ko','Markersize',5)

subplot(5,2,7)
plot(CFDataWT(1:20,4),'-k')
title('WT 50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50_WT(1:20),'ko','Markersize',5)

subplot(5,2,9)
semilogx(rec,CFDataWT(21:30,4),'-k')
title('WT 50 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_50_WT(21:30),'ko','Markersize',5)

subplot(5,2,2)
plot(CFDataKO(1:20,1),'-k')
title('KO 1 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_1_KO,'ko','Markersize',5)

subplot(5,2,4)
plot(CFDataKO(1:20,2),'-k')
title('KO 10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10_KO,'ko','Markersize',5)

subplot(5,2,6)
plot(CFDataKO(1:20,3),'-k')
title('KO 20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20_KO,'ko','Markersize',5)

subplot(5,2,8)
plot(CFDataKO(1:20,4),'-k')
title('KO 50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50_KO(1:20),'ko','Markersize',5)

subplot(5,2,10)
semilogx(rec,CFDataKO(21:30,4),'-k')
title('KO 50 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_50_KO(21:30),'ko','Markersize',5)

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, SS, reserve_size, k_refill, Syt3)

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
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt3(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [~,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill, Syt3(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,reserve_size,k_refill,Syt3)

multFactor = 1+Syt3;
dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*multFactor + state(2)*k_undocking*multFactor;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;

end

function dydt = dState(t,state,k_docking,k_undocking,reserve_size,k_refill,Syt3)

multFactor = 1+Syt3(round(t/.1)+1);
dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*multFactor + state(2)*k_undocking*multFactor;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;
end