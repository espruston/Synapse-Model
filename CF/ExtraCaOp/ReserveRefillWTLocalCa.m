CFDataWT = load('..\CFDataWT.mat').CFDataWT;
CF1HzCa = load('..\CF1HzCa.mat').Ca_sim;
CF10HzCa = load('..\CF10HzCa.mat').Ca_sim;
CF20HzCa = load('..\CF20HzCa.mat').Ca_sim;
CF50HzCa = load('..\CF50HzCa.mat').Ca_sim;

%x = [.95,0.004,0.002,2.5,0.00088,1e6];
%x = [0.7304,0.0024994,0.0005734,8.479,0.0001,100000.0001]; %best fit
%1/17/2021 err=0.764
%x = [0.70303,0.0033015,0.00010035,6.7483,0.0001,0]; %KO best fit 1/13/2021
%err = 0.69845
x = [0.70303,0.0033015,0.00010035,6.7483,0.0001,35]; 

p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);
k_refill = x(5);

C_3 = x(6); %multiplicative factor of syt3 effect

k_on_3 = 0.3; %ms^-1 Hui
k_off_3 = 0.1929; %ms^-1  Hui
Ca_rest = 50e-9;

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
state_0 = [1; 0; reserve_size; 0]; %[empty pool; bound pool; reserve pool; Syt3]

[t0,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, Ca_rest, reserve_size, k_refill, C_3, k_on_3, k_off_3), [0 t_SS], state_0);

SS = state(end,:);

[ts_1, state_1, Fused_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, CF1HzCa, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_1 = Fused_1/Fused_1(1);

[ts_10, state_10, Fused_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, CF10HzCa, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_10 = Fused_10/Fused_10(1);

[ts_20, state_20, Fused_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, CF20HzCa, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_20 = Fused_20/Fused_20(1);

[ts_50, state_50, Fused_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, CF50HzCa, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_50 = Fused_50/Fused_50(1);

%err_3KO = sqrt(sum((hz_1 - CFData3KO(1:20,1)).^2 + (hz_10 - CFData3KO(1:20,2)).^2 + (hz_20 - CFData3KO(1:20,3)).^2 + (hz_50(1:20) - CFData3KO(1:20,4)).^2) + 10*sum((hz_50(21:30) - CFData3KO(21:30,4)).^2));
err_WT = sqrt(2*sum((hz_1(1:5) - CFDataWT(1:5,1)).^2 + (hz_10(1:5) - CFDataWT(1:5,2)).^2 + (hz_20(1:5) - CFDataWT(1:5,3)).^2 + (hz_50(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_1(6:20) - CFDataWT(6:20,1)).^2 + (hz_10(6:20) - CFDataWT(6:20,2)).^2 + (hz_20(6:20) - CFDataWT(6:20,3)).^2 + (hz_50(6:20) - CFDataWT(6:20,4)).^2) + 10*sum((hz_50(21:30) - CFDataWT(21:30,4)).^2));

err = err_WT;

rec = [50 100 200 350 500 750 1000 2000 5000 10000];

figure('Name','WT CF Simulations','NumberTitle','off')
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

figure('Name','WT CF Simulations','NumberTitle','off')
subplot(2,2,1)
plot(ts_1,state_1(:,1),'color',[255, 165, 0]/255)
title('1 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_1,state_1(:,2),'-b')
plot(ts_1,state_1(:,3),'-k')

ax = gca;
ax.XRuler.Exponent = 0;

legend({'Empty Sites','Filled Sites','Reserve Vesicles'},'Location','Best')

subplot(2,2,2)
plot(ts_10,state_10(:,1),'color',[255, 165, 0]/255)
title('10 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_10,state_10(:,2),'-b')
plot(ts_10,state_10(:,3),'-k')

subplot(2,2,3)
plot(ts_20,state_20(:,1),'color',[255, 165, 0]/255)
title('20 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_20,state_20(:,2),'-b')
plot(ts_20,state_20(:,3),'-k')

subplot(2,2,4)
plot(ts_50,state_50(:,1),'color',[255, 165, 0]/255)
title('50 Hz')
xlabel('time (ms)')
ylabel('Value')
set(gca,'ylim',[0 reserve_size])
hold on
plot(ts_50,state_50(:,2),'-b')
plot(ts_50,state_50(:,3),'-k')

%plot Ca and Syts
figure('Name','WT CF Simulations','NumberTitle','off')
subplot(4,2,1)
semilogy(ts_1,CF1HzCa(round(ts_1/.1)+1))
title('1 Hz Ca')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(4,2,2)
plot(ts_1,state_1(:,4))
title('1 Hz Syt3')
xlabel('time (ms)')
ylabel('Bound Syt3')

subplot(4,2,3)
semilogy(ts_10,CF10HzCa(round(ts_10/.1)+1))
title('10 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(4,2,4)
plot(ts_10,state_10(:,4))
title('10 Hz')
xlabel('time (ms)')
ylabel('Bound Syt3')

subplot(4,2,5)
semilogy(ts_20,CF20HzCa(round(ts_20/.1)+1))
title('20 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(4,2,6)
plot(ts_20,state_20(:,4))
title('20 Hz')
xlabel('time (ms)')
ylabel('Bound Syt3')

subplot(4,2,7)
semilogy(ts_50,CF50HzCa(1,round(ts_50/.1)+1))
title('50 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (M)')

subplot(4,2,8)
plot(ts_50,state_50(:,4))
title('50 Hz')
xlabel('time (ms)')
ylabel('Bound Syt3')

disp(err)

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, Ca, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release, -pre_stim(2)*p_release, 0, 0];
        Fused(i) = pre_stim(2)*p_release;  
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, Ca(1,:), reserve_size, k_refill, C_3, k_on_3, k_off_3), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, Ca(i,:), reserve_size, k_refill, C_3, k_on_3, k_off_3), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,Ca_rest,reserve_size,k_refill,C_3,k_on_3,k_off_3)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking + state(2)*k_undocking*C_3*(1-state(4));
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;
dydt(4,1) = (1-state(4))*k_on_3*(Ca_rest^3.6/(Ca_rest^3.6+7e-6^3.6)) - state(4)*k_off_3;

end

function dydt = dState(t,state,k_docking,k_undocking,Ca,reserve_size,k_refill,C_3,k_on_3,k_off_3)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking + state(2)*k_undocking*C_3*(1-state(4));
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;
dydt(4,1) = (1-state(4))*k_on_3*(Ca(round(t/.1)+1)^3.6/(Ca(round(t/.1)+1)^3.6+7e-6^3.6)) - state(4)*k_off_3;

end