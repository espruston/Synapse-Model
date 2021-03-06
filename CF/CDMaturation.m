%x = [0.57676,0.063374,0.0043336,0.0030516,0.004604,0.00118,3.8875];
%x = [0.95,0,0.0019237,0.0006432,0.0037323,0.0077772,10];
x = [0.73335,0.11204,0.0038305,0.0037433,0.0064804,0.001609,9.0626];

p_mature = x(1); 
p_immature = x(2); 
k_docking = x(3); 
k_undocking = x(4); 
k_maturation = x(5); 
k_dematuration = x(6); 
C_3 = x(7);

C_7 = 0;

%Load CF data from 200916
%CFHz = [1 10 20 50];
%CFRecovery = [50 100 200 350 500 750 1000 2000 5000 10000];
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

state_0 = [1; 0; 0;]; %[empty; immature; mature]

[~,state] = ode113(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation*(1+CF1HzSyt3(1)*C_3),k_dematuration), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_im_1, Fused_m_1] = stim_sim(stimulus_times_1, max_time_1, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF1HzSyt3, C_7, CF1HzSyt7, SS);
Fused_1 = Fused_im_1 + Fused_m_1;
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_im_10, Fused_m_10] = stim_sim(stimulus_times_10, max_time_10, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF10HzSyt3, C_7, CF10HzSyt7, SS);
Fused_10 = Fused_im_10 + Fused_m_10;
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_im_20, Fused_m_20] = stim_sim(stimulus_times_20, max_time_20, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF20HzSyt3, C_7, CF20HzSyt7, SS);
Fused_20 = Fused_im_20 + Fused_m_20;
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_im_50, Fused_m_50] = stim_sim(stimulus_times_50, max_time_50, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF50HzSyt3, C_7, CF50HzSyt7, SS);
Fused_50 = Fused_im_50 + Fused_m_50;
hz_50 = Fused_50/Fused_50(1);

err_WT = sqrt(sum((hz_1 - CFDataWT(1:20,1)).^2 + (hz_10 - CFDataWT(1:20,2)).^2 + (hz_20(1:20) - CFDataWT(1:20,3)).^2 + (hz_50 - CFDataWT(1:20,4)).^2) + sum(10*(hz_20(21:30) - CFDataWT(21:30,3)).^2));

Fused_1_norm = Fused_im_1(1) + Fused_m_1(1);
Fused_10_norm = Fused_im_10(1) + Fused_m_10(1);
Fused_20_norm = Fused_im_20(1) + Fused_m_20(1);
Fused_50_norm = Fused_im_50(1) + Fused_m_50(1);

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
plot(Fused_im_1/Fused_1_norm,'rv')
plot(Fused_m_1/Fused_1_norm,'g^')
%legend(labels(1),'1 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')


subplot(3,2,2)
plot(CFDataWT(1:20,2),'-k')
title('10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10,'ko','Markersize',5)
plot(Fused_im_10/Fused_10_norm,'rv')
plot(Fused_m_10/Fused_10_norm,'g^')

subplot(3,2,3)
plot(CFDataWT(1:20,3),'-k')
title('20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20(1:20),'ko','Markersize',5)
plot(Fused_im_20(1:20)/Fused_20_norm,'rv')
plot(Fused_m_20(1:20)/Fused_20_norm,'g^')

subplot(3,2,4)
plot(CFDataWT(1:20,4),'-k')
title('50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50,'ko','Markersize',5)
plot(Fused_im_50/Fused_50_norm,'rv')
plot(Fused_m_50/Fused_50_norm,'g^')

subplot(3,2,5)
semilogx(rec,CFDataWT(21:30,3),'-k')
title('20 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_20(21:30),'ko','Markersize',5)
semilogx(rec,Fused_im_20(21:30)/Fused_20_norm,'rv')
semilogx(rec,Fused_m_20(21:30)/Fused_20_norm,'g^')

%Syt3 KO
C_3 = 0;
[t0,state] = ode113(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation*(1+CF1HzSyt3(1)*C_3),k_dematuration), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_im_1, Fused_m_1] = stim_sim(stimulus_times_1, max_time_1, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF1HzSyt3, C_7, CF1HzSyt7, SS);
Fused_1 = Fused_im_1 + Fused_m_1;
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_im_10, Fused_m_10] = stim_sim(stimulus_times_10, max_time_10, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF10HzSyt3, C_7, CF10HzSyt7, SS);
Fused_10 = Fused_im_10 + Fused_m_10;
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_im_20, Fused_m_20] = stim_sim(stimulus_times_20, max_time_20, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF20HzSyt3, C_7, CF20HzSyt7, SS);
Fused_20 = Fused_im_20 + Fused_m_20;
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_im_50, Fused_m_50] = stim_sim(stimulus_times_50, max_time_50, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, CF50HzSyt3, C_7, CF50HzSyt7, SS);
Fused_50 = Fused_im_50 + Fused_m_50;
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
plot(Fused_im_1/Fused_1_norm,'rv')
plot(Fused_m_1/Fused_1_norm,'g^')
%legend(labels(1),'1 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')


subplot(3,2,2)
plot(CFData3KO(1:20,2),'-k')
title('10 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1])
hold on
plot(hz_10,'ko','Markersize',5)
plot(Fused_im_10/Fused_10_norm,'rv')
plot(Fused_m_10/Fused_10_norm,'g^')

subplot(3,2,3)
plot(CFData3KO(1:20,3),'-k')
title('20 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_20(1:20),'ko','Markersize',5)
plot(Fused_im_20(1:20)/Fused_20_norm,'rv')
plot(Fused_m_20(1:20)/Fused_20_norm,'g^')

subplot(3,2,4)
plot(CFData3KO(1:20,4),'-k')
title('50 Hz')
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 20])
set(gca,'ylim',[0 1.2])
hold on
plot(hz_50,'ko','Markersize',5)
plot(Fused_im_50/Fused_50_norm,'rv')
plot(Fused_m_50/Fused_50_norm,'g^')

subplot(3,2,5)
semilogx(rec,CFData3KO(21:30,3),'-k')
title('20 Hz Recovery')
xlabel('t (ms)')
ylabel('Peak EPSC')
set(gca,'xlim',[50 10000])
set(gca,'ylim',[0 1.2])
hold on
semilogx(rec,hz_20(21:30),'ko','Markersize',5)
semilogx(rec,Fused_im_20(21:30)/Fused_20_norm,'rv')
semilogx(rec,Fused_m_20(21:30)/Fused_20_norm,'g^')

err = (err_WT + err_3KO)/2;

disp(err)

function [ts, state, Fused_im, Fused_m] = stim_sim(stimulus_times, max_time, p_mature, p_immature, k_docking, k_undocking, k_maturation, k_dematuration, C_3, Syt3, C_7, Syt7, SS)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused_im = zeros(length(stimulus_times),1);
    Fused_m = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_immature+pre_stim(3)*p_mature*(1+C_7*Syt7(1,stimulus_times(i)/0.01+1)), -pre_stim(2)*p_immature, -pre_stim(3)*p_mature*(1+C_7*Syt7(1,stimulus_times(i)/0.01+1))];
        Fused_im(i) = pre_stim(2)*p_immature;
        Fused_m(i) = pre_stim(3)*p_mature*(1+C_7*Syt7(1,stimulus_times(i)/0.01+1));
        [t,out] = ode113(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,C_3,Syt3(1,:)), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
    if stimulus_times == linspace(0,50*19,20)
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_im_rec = zeros(length(Recovery),1);
        Fused_m_rec = zeros(length(Recovery),1);
        
        for i = 1:length(Recovery)
            
            [t,out] = ode113(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,C_3,Syt3(i,:)), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_im_rec(i) = pre_stim(2)*p_immature;
            Fused_m_rec(i) = pre_stim(3)*p_mature*(1+C_7*Syt7(stimulus_times(i)/0.01+1));
            
        end
        Fused_im = [Fused_im; Fused_im_rec];
        Fused_m = [Fused_m; Fused_m_rec];
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end
    
function dydt = dSS(~,state,k_docking,k_undocking,k_maturation,k_dematuration)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,C_3,Syt3)
    
    k_maturation = k_maturation*(1+Syt3(round(t/.01)+1)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end