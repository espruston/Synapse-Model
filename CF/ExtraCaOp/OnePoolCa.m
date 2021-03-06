x = [0.85,0.001,0.0001,50,0,7e6];
p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);

C_3 = x(5);
C_Ca = x(6);

Ca_rest = 50e-9;

if exist('loaded','var') == 0

    %Load CF data from 200916
    CFDataWT = load('..\CFDataWT.mat').CFDataWT;
    CFData3KO = load('..\CFData3KO.mat').CFData3KO;

    %load CF syt sims
    CF1HzSyt3 = load('..\CF1HzSyt3.mat').Syt3;
    CF10HzSyt3 = load('..\CF10HzSyt3.mat').Syt3;
    CF20HzSyt3 = load('..\CF20HzSyt3.mat').Syt3;
    CF50HzSyt3 = load('..\CF50HzSyt3.mat').Syt3;

    %load Ca Signals
    CF1HzCa = load('..\CF1HzCa.mat').Ca_sim;
    CF10HzCa = load('..\CF10HzCa.mat').Ca_sim;
    CF20HzCa = load('..\CF20HzCa.mat').Ca_sim;
    CF50HzCa = load('..\CF50HzCa.mat').Ca_sim;
    
    loaded = 1;
    
end

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

%Check WT
[~,state] = ode45(@(t,state) dSS(t, state, k_docking, k_undocking, CF1HzSyt3(1)*C_3, Ca_rest*C_Ca, reserve_size), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, CF1HzSyt3*C_3, CF1HzCa*C_Ca, SS, reserve_size);
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, CF10HzSyt3*C_3, CF10HzCa*C_Ca, SS, reserve_size);
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, CF20HzSyt3*C_3, CF20HzCa*C_Ca, SS, reserve_size);
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, CF50HzSyt3*C_3, CF50HzCa*C_Ca, SS, reserve_size);
hz_50 = Fused_50/Fused_50(1);

err_WT = sqrt(sum((hz_1 - CFDataWT(1:20,1)).^2 + (hz_10 - CFDataWT(1:20,2)).^2 + (hz_20 - CFDataWT(1:20,3)).^2 + (hz_50(1:20) - CFDataWT(1:20,4)).^2) + sum(10*(hz_50(21:30) - CFDataWT(21:30,4)).^2));

rec = [50 100 200 350 500 750 1000 2000 5000 10000];

figure('Name','WT Simulations','NumberTitle','off')
subplot(3,2,1)
plot(CFDataWT(1:20,1),'-k')
title('1 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_1,'ko','Markersize',5)
%legend({'Data','Simulation'},'Position',[0.14, 0.75, .1, .1])

subplot(3,2,2)
plot(CFDataWT(1:20,2),'-k')
title('10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10,'ko','Markersize',5)

subplot(3,2,3)
plot(CFDataWT(1:20,3),'-k')
title('20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20,'ko','Markersize',5)

subplot(3,2,4)
plot(CFDataWT(1:20,4),'-k')
title('50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50(1:20),'ko','Markersize',5)

subplot(3,2,5)
semilogx(rec,CFDataWT(21:30,4),'-k')
title('50 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_50(21:30),'ko','Markersize',5)

%Syt3 KO
C_3 = 0;
[t0,state] = ode45(@(t,state) dSS(t, state, k_docking, k_undocking, C_3, Ca_rest*C_Ca, reserve_size), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, CF1HzSyt3*C_3, CF1HzCa*C_Ca, SS, reserve_size);
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, CF10HzSyt3*C_3, CF10HzCa*C_Ca, SS, reserve_size);
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, CF20HzSyt3*C_3, CF20HzCa*C_Ca, SS, reserve_size);
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, CF50HzSyt3*C_3, CF50HzCa*C_Ca, SS, reserve_size);
hz_50 = Fused_50/Fused_50(1);

err_3KO = sqrt(sum((hz_1 - CFData3KO(1:20,1)).^2 + (hz_10 - CFData3KO(1:20,2)).^2 + (hz_20 - CFData3KO(1:20,3)).^2 + (hz_50(1:20) - CFData3KO(1:20,4)).^2) + sum(10*(hz_50(21:30) - CFData3KO(21:30,4)).^2));

err = (err_WT + err_3KO)/2;

figure('Name','Syt3 KO Simulations','NumberTitle','off')
subplot(3,2,1)
plot(CFData3KO(1:20,1),'-k')
title('1 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_1,'ko','Markersize',5)
legend({'Data','Simulation'},'Position',[0.14, 0.75, .1, .1])

subplot(3,2,2)
plot(CFData3KO(1:20,2),'-k')
title('10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10,'ko','Markersize',5)

subplot(3,2,3)
plot(CFData3KO(1:20,3),'-k')
title('20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20,'ko','Markersize',5)

subplot(3,2,4)
plot(CFData3KO(1:20,4),'-k')
title('50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50(1:20),'ko','Markersize',5)

subplot(3,2,5)
semilogx(rec,CFData3KO(21:30,4),'-k')
title('50 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_50(21:30),'ko','Markersize',5)

disp(err)

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, Syt3, Ca, SS, reserve_size)

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
        [t,out] = ode45(@(t,state) dState(t, state, k_docking, k_undocking, Syt3(1,:), Ca(1,:), reserve_size), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20)
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [t,out] = ode45(@(t,state) dState(t, state, k_docking, k_undocking, Syt3(i,:), Ca(i,:), reserve_size), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,Syt3,Ca,reserve_size)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+Syt3) - state(1)*(state(3)/reserve_size)*k_docking*(1+Ca) + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1);

end

function dydt = dState(t,state,k_docking,k_undocking,Syt3,Ca,reserve_size)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+Syt3(round(t/.01)+1)) - state(1)*(state(3)/reserve_size)*k_docking*(1+Ca(round(t/.01)+1)) + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1);

end

%interdependent ca mechanism & syt3

% function dydt = dSS(~,state,k_docking,k_undocking,Syt3,Ca,reserve_size)
% 
% dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+Syt3)*(1+Ca) + state(2)*k_undocking;
% dydt(2,1) = -dydt(1,1);
% dydt(3,1) = dydt(1,1);
% 
% end
% 
% function dydt = dState(t,state,k_docking,k_undocking,Syt3,Ca,reserve_size)
% 
% dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+Syt3(round(t/.01)+1))*(1+Ca(round(t/.01)+1)) + state(2)*k_undocking;
% dydt(2,1) = -dydt(1,1);
% dydt(3,1) = dydt(1,1);
% 
% end
