%x = [0.68146,0.00035437,0.0027224,0.0033204,0.0079877,0.49609,10]; %best fit 11/16/2020 err = 1.0954 
x = [0.77081,0.0013941,0.0013708,0.0021637,0.0032797,0.59978,10]; %best fit 11/23/2020 err=1.0161

p_release = x(1); 
k_docking_1 = x(2); 
k_undocking_1 = x(3); 
k_docking_2 = x(4); 
k_undocking_2 = x(5); 
size_rel = x(6);
C_3 = x(7);

C_7 = 0;

%Load CF data from 200916
CFHz = [1 10 20 50];
CFRecovery = [50 100 200 350 500 750 1000 2000 5000 10000];
CFDataWT = matfile('CFDataWT.mat').CFDataWT;
CFData3KO = matfile('CFData3KO.mat').CFData3KO;
CFData7KO = matfile('CFData7KO.mat').CFData7KO;

%load CF syt sims
CF1HzSyt3 = matfile('CF1HzSyt3.mat').Syt3;
CF10HzSyt3 = matfile('CF10HzSyt3.mat').Syt3;
CF20HzSyt3 = matfile('CF20HzSyt3.mat').Syt3;
CF50HzSyt3 = matfile('CF50HzSyt3.mat').Syt3;
CF1HzSyt7 = matfile('CF1HzSyt7.mat').Syt7;
CF10HzSyt7 = matfile('CF10HzSyt7.mat').Syt7;
CF20HzSyt7 = matfile('CF20HzSyt7.mat').Syt7;
CF50HzSyt7 = matfile('CF50HzSyt7.mat').Syt7;

%Define train stims
stimulus_times_1 = linspace(0,1000*19,20); %100 stims 1hz
max_time_1 = stimulus_times_1(end) + stimulus_times_1(2)*3;

stimulus_times_10 = linspace(0,100*19,20); %100 stims 10hz
max_time_10 = stimulus_times_10(end) + stimulus_times_10(2)*3;

stimulus_times_20 = linspace(0,50*19,20);
max_time_20 = stimulus_times_20(end) + stimulus_times_20(2)*3;
%stimulus_times_20 = [stimulus_times_20 stimulus_times_20(end)+[50 100 200 350 500 750 1000 2000 5000 10000]];

stimulus_times_50 = linspace(0,20*19,20);
max_time_50 = stimulus_times_50(end) + stimulus_times_50(2)*3;

%Check WT
t_SS = 10000; %ms

state_0 = [size_rel; 1-size_rel; 0; 0]; %[empty pool 1; empty pool 2; pool 1; pool 2 (syt3 pool)]

[~,state] = ode113(@(t,state) dSS(t,state,k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF1HzSyt3(1)*C_3), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_2_1, Fused_1_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF1HzSyt3*C_3, C_7, CF1HzSyt7, SS);
Fused_1 = Fused_2_1 + Fused_1_1;
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_2_10, Fused_1_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF10HzSyt3*C_3, C_7, CF10HzSyt7, SS);
Fused_10 = Fused_2_10 + Fused_1_10;
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_2_20, Fused_1_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF20HzSyt3*C_3, C_7, CF20HzSyt7, SS);
Fused_20 = Fused_2_20 + Fused_1_20;
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_2_50, Fused_1_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF50HzSyt3*C_3, C_7, CF50HzSyt7, SS);
Fused_50 = Fused_2_50 + Fused_1_50;
hz_50 = Fused_50/Fused_50(1);

err_WT = sqrt(sum((hz_1 - CFDataWT(1:20,1)).^2 + (hz_10 - CFDataWT(1:20,2)).^2 + (hz_20(1:20) - CFDataWT(1:20,3)).^2 + (hz_50 - CFDataWT(1:20,4)).^2) + sum(10*(hz_20(21:30) - CFDataWT(21:30,3)).^2));

Fused_1_norm = Fused_2_1(1) + Fused_1_1(1);
Fused_10_norm = Fused_2_10(1) + Fused_1_10(1);
Fused_20_norm = Fused_2_20(1) + Fused_1_20(1);
Fused_50_norm = Fused_2_50(1) + Fused_1_50(1);

labels = ["1 Hz data", "10 Hz data", "20 Hz data", "50 Hz data"];
rec = [50 100 200 350 500 750 1000 2000 5000 10000];

figure 
subplot(3,2,1)
plot(CFDataWT(1:20,1),'-k')
title('1 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_1,'ko','Markersize',5)
plot(Fused_2_1/Fused_1_norm,'rv')
plot(Fused_1_1/Fused_1_norm,'g^')
%legend(labels(1),'1 Hz Simulation','Syt3 Pool Fusion','Syt3 Ind. Pool Fusion')


subplot(3,2,2)
plot(CFDataWT(1:20,2),'-k')
title('10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10,'ko','Markersize',5)
plot(Fused_2_10/Fused_10_norm,'rv')
plot(Fused_1_10/Fused_10_norm,'g^')

subplot(3,2,3)
plot(CFDataWT(1:20,3),'-k')
title('20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20(1:20),'ko','Markersize',5)
plot(Fused_2_20(1:20)/Fused_20_norm,'rv')
plot(Fused_1_20(1:20)/Fused_20_norm,'g^')

subplot(3,2,4)
plot(CFDataWT(1:20,4),'-k')
title('50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50,'ko','Markersize',5)
plot(Fused_2_50/Fused_50_norm,'rv')
plot(Fused_1_50/Fused_50_norm,'g^')

subplot(3,2,5)
semilogx(rec,CFDataWT(21:30,3),'-k')
title('20 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_20(21:30),'ko','Markersize',5)
semilogx(rec,Fused_2_20(21:30)/Fused_20_norm,'rv')
semilogx(rec,Fused_1_20(21:30)/Fused_20_norm,'g^')

%Syt3 KO
C_3 = 0;
[t0,state] = ode113(@(t,state) dSS(t,state,k_docking_1, k_undocking_1, k_docking_2, k_undocking_2,C_3), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_2_1, Fused_1_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF1HzSyt3*C_3, C_7, CF1HzSyt7, SS);
Fused_1 = Fused_2_1 + Fused_1_1;
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_2_10, Fused_1_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF10HzSyt3*C_3, C_7, CF10HzSyt7, SS);
Fused_10 = Fused_2_10 + Fused_1_10;
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_2_20, Fused_1_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF20HzSyt3*C_3, C_7, CF20HzSyt7, SS);
Fused_20 = Fused_2_20 + Fused_1_20;
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_2_50, Fused_1_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, CF50HzSyt3*C_3, C_7, CF50HzSyt7, SS);
Fused_50 = Fused_2_50 + Fused_1_50;
hz_50 = Fused_50/Fused_50(1);

err_3KO = sqrt(sum((hz_1 - CFData3KO(1:20,1)).^2 + (hz_10 - CFData3KO(1:20,2)).^2 + (hz_20(1:20) - CFData3KO(1:20,3)).^2 + (hz_50 - CFData3KO(1:20,4)).^2) + sum(10*(hz_20(21:30) - CFData3KO(21:30,3)).^2));

figure
subplot(3,2,1)
plot(CFData3KO(1:20,1),'-k')
title('1 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_1,'ko','Markersize',5)
plot(Fused_2_1/Fused_1_norm,'rv')
plot(Fused_1_1/Fused_1_norm,'g^')
%legend(labels(1),'1 Hz Simulation','Syt3 Pool Fusion','Syt3 Ind. Pool Fusion')


subplot(3,2,2)
plot(CFData3KO(1:20,2),'-k')
title('10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10,'ko','Markersize',5)
plot(Fused_2_10/Fused_10_norm,'rv')
plot(Fused_1_10/Fused_10_norm,'g^')

subplot(3,2,3)
plot(CFData3KO(1:20,3),'-k')
title('20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20(1:20),'ko','Markersize',5)
plot(Fused_2_20(1:20)/Fused_20_norm,'rv')
plot(Fused_1_20(1:20)/Fused_20_norm,'g^')

subplot(3,2,4)
plot(CFData3KO(1:20,4),'-k')
title('50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50,'ko','Markersize',5)
plot(Fused_2_50/Fused_50_norm,'rv')
plot(Fused_1_50/Fused_50_norm,'g^')

subplot(3,2,5)
semilogx(rec,CFData3KO(21:30,3),'-k')
title('20 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_20(21:30),'ko','Markersize',5)
semilogx(rec,Fused_2_20(21:30)/Fused_20_norm,'rv')
semilogx(rec,Fused_1_20(21:30)/Fused_20_norm,'g^')

err = (err_WT + err_3KO)/2;

disp(err)

function [ts, state, Fused_2, Fused_1] = stim_sim(stimulus_times, max_time, p_release, k_docking_1, k_undocking_1, k_docking_2, k_undocking_2, Syt3, C_7, Syt7, SS)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused_1 = zeros(length(stimulus_times),1);
    Fused_2 = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(3)*p_release, pre_stim(4)*p_release, -pre_stim(3)*p_release, -pre_stim(4)*p_release];
        Fused_1(i) = pre_stim(3)*p_release;
        Fused_2(i) = pre_stim(4)*p_release;   
        [t,out] = ode113(@(t,state) dState(t,state,k_docking_1, k_undocking_1, k_docking_2, k_undocking_2,Syt3(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
    if stimulus_times == linspace(0,50*19,20)
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_1_rec = zeros(length(Recovery),1);
        Fused_2_rec = zeros(length(Recovery),1);    
        
        for i = 1:length(Recovery)
            
            [t,out] = ode113(@(t,state) dState(t,state,k_docking_1, k_undocking_1, k_docking_2, k_undocking_2,Syt3(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_1_rec(i) = pre_stim(3)*p_release;
            Fused_2_rec(i) = pre_stim(4)*p_release;
       
        end
        Fused_1 = [Fused_1; Fused_1_rec];
        Fused_2 = [Fused_2; Fused_2_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end
    
function dydt = dSS(~,state,k_docking_1,k_undocking_1,k_docking_2,k_undocking_2,Syt3)
    
    dydt(1,1) = -state(1)*k_docking_1 + state(3)*k_undocking_1;
    dydt(2,1) = -state(2)*k_docking_2*Syt3 + state(4)*k_undocking_2;
    dydt(3,1) = -dydt(1,1);
    dydt(4,1) = -dydt(2,1);

end

function dydt = dState(t,state,k_docking_1,k_undocking_1,k_docking_2,k_undocking_2,Syt3)
    
    dydt(1,1) = -state(1)*k_docking_1 + state(3)*k_undocking_1;
    dydt(2,1) = -state(2)*k_docking_2*Syt3(round(t/.01)+1) + state(4)*k_undocking_2;
    dydt(3,1) = -dydt(1,1);
    dydt(4,1) = -dydt(2,1);

end