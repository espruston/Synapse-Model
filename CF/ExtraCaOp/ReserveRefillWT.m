CFDataWT = load('..\CFDataWT.mat').CFDataWT;

%x = [.95,0.004,0.002,2.5,0.00088,1e6];
%x = [0.7304,0.0024994,0.0005734,8.479,0.0001,100000.0001]; %best fit
%1/17/2021 err=0.764
x = [0.70303,0.0033015,0.00010035,6.7483,0.0001,0]; %best fit 1/13/2021
%err = 0.69845
 
p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);
k_refill = x(5);

C_3 = x(6); %multiplicative factor of syt3 effect

% k_on_3 = 3e5; %M^-1ms^-1 Hui
% k_off_3 = 0.2; %ms^-1  Hui
k_on_3 = .157; %M^-1ms^-1 Hui
k_off_3 = 0.011; %ms^-1  Hui
Ca_rest = 50e-9;
Ca_residual = 250e-9;
k_Ca_decay = 1/48;

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
state_0 = [1; 0; reserve_size; Ca_rest; 0]; %[empty pool; bound pool; reserve pool; Ca[M]; Syt3]

%Syt3 KO
[t0,state] = ode15s(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill, C_3, k_on_3, k_off_3), [0 t_SS], state_0);

SS = state(end,:);

[ts_1, state_1, Fused_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, Ca_rest, Ca_residual, k_Ca_decay, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_1 = Fused_1/Fused_1(1);

[ts_10, state_10, Fused_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, Ca_rest, Ca_residual, k_Ca_decay, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_10 = Fused_10/Fused_10(1);

[ts_20, state_20, Fused_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, Ca_rest, Ca_residual, k_Ca_decay, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_20 = Fused_20/Fused_20(1);

[ts_50, state_50, Fused_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, Ca_rest, Ca_residual, k_Ca_decay, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3);
hz_50 = Fused_50/Fused_50(1);

%err_3KO = sqrt(sum((hz_1 - CFData3KO(1:20,1)).^2 + (hz_10 - CFData3KO(1:20,2)).^2 + (hz_20 - CFData3KO(1:20,3)).^2 + (hz_50(1:20) - CFData3KO(1:20,4)).^2) + 10*sum((hz_50(21:30) - CFData3KO(21:30,4)).^2));
err_3KO = sqrt(2*sum((hz_1(1:5) - CFDataWT(1:5,1)).^2 + (hz_10(1:5) - CFDataWT(1:5,2)).^2 + (hz_20(1:5) - CFDataWT(1:5,3)).^2 + (hz_50(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_1(6:20) - CFDataWT(6:20,1)).^2 + (hz_10(6:20) - CFDataWT(6:20,2)).^2 + (hz_20(6:20) - CFDataWT(6:20,3)).^2 + (hz_50(6:20) - CFDataWT(6:20,4)).^2) + 10*sum((hz_50(21:30) - CFDataWT(21:30,4)).^2));

err = err_3KO;

rec = [50 100 200 350 500 750 1000 2000 5000 10000];

figure('Name','Syt3 KO Simulations','NumberTitle','off')
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

figure('Name','Syt3 KO Simulations','NumberTitle','off')
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
figure
subplot(4,2,1)
plot(ts_1,state_1(:,4)*1e9)
title('1 Hz Ca')
xlabel('time (ms)')
ylabel('Ca Conc. (nM)')
set(gca,'ylim',[0 300])

ax = gca;
ax.XRuler.Exponent = 0;


subplot(4,2,2)
%plot(ts_10,state_1(:,5))
plot(ts_1,state_1(:,4).^3.6./(state_1(:,4).^3.6+7e-6^3.6))
title('1 Hz Syt3')
xlabel('time (ms)')
ylabel('Bound Syt3')

subplot(4,2,3)
plot(ts_10,state_10(:,4)*1e9)
title('10 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (nM)')
set(gca,'ylim',[0 300])

subplot(4,2,4)
plot(ts_10,state_10(:,5))
title('10 Hz')
xlabel('time (ms)')
ylabel('Bound Syt3')

subplot(4,2,5)
plot(ts_20,state_20(:,4)*1e9)
title('20 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (nM)')
set(gca,'ylim',[0 500])

subplot(4,2,6)
plot(ts_20,state_20(:,5))
title('20 Hz')
xlabel('time (ms)')
ylabel('Bound Syt3')

subplot(4,2,7)
plot(ts_50,state_50(:,4)*1e9)
title('50 Hz')
xlabel('time (ms)')
ylabel('Ca Conc. (nM)')
set(gca,'ylim',[0 300])

subplot(4,2,8)
plot(ts_50,state_50(:,5))
title('50 Hz')
xlabel('time (ms)')
ylabel('Bound Syt3')

disp(err)

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, Ca_rest, Ca_residual, k_Ca_decay, SS, reserve_size, k_refill, C_3, k_on_3, k_off_3)

    delta_t = 1e-2; %ms

    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    state = SS;
    
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release, -pre_stim(2)*p_release, 0, (Ca_residual-pre_stim(4)), 0];
        Fused(i) = pre_stim(2)*p_release;  
        [t,out] = ode15s(@(t,state) dState(t, state, k_docking, k_undocking, Ca_rest, k_Ca_decay, reserve_size, k_refill, C_3, k_on_3, k_off_3), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [t,out] = ode15s(@(t,state) dState(t, state, k_docking, k_undocking, Ca_rest, k_Ca_decay, reserve_size, k_refill, C_3, k_on_3, k_off_3), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,reserve_size,k_refill,C_3,k_on_3,k_off_3)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+C_3*state(5)) + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;
dydt(4,1) = 0;
dydt(5,1) = (1-state(5))*k_on_3*(state(4)^3.6/(state(4)^3.6+7e-6^3.6)) - state(5)*k_off_3;

end

function dydt = dState(t,state,k_docking,k_undocking,Ca_rest,k_Ca_decay,reserve_size,k_refill,C_3,k_on_3,k_off_3)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking*(1+C_3*state(5)) + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;
dydt(4,1) = -k_Ca_decay*(state(4)-Ca_rest);
dydt(5,1) = (1-state(5))*k_on_3*(state(4)^3.6/(state(4)^3.6+7e-6^3.6)) - state(5)*k_off_3;

end