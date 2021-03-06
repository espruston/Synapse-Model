function err = ReserveRefillFunc(x)

CFData3KO = load('..\CFData3KO.mat').CFData3KO;
CFDataWT = load('..\CFDataWT.mat').CFDataWT;
CF1HzCa = load('..\CF1HzCa.mat').Ca_sim;
CF10HzCa = load('..\CF10HzCa.mat').Ca_sim;
CF20HzCa = load('..\CF20HzCa.mat').Ca_sim;
CF50HzCa = load('..\CF50HzCa.mat').Ca_sim;

p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);
k_refill = x(5);

C_3 = x(6);

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
state_0 = [1; 0; reserve_size]; %[empty pool; bound pool; reserve pool]

%Syt3 KO
[t0,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, Ca_rest, reserve_size, k_refill), [0 t_SS], state_0);

SS = state(end,:);

[~, ~, Fused_1] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, CF1HzCa, SS, reserve_size, k_refill);
hz_1 = Fused_1/Fused_1(1);

[~, ~, Fused_10] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, CF10HzCa, SS, reserve_size, k_refill);
hz_10 = Fused_10/Fused_10(1);

[~, ~, Fused_20] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, CF20HzCa*C_Ca, SS, reserve_size, k_refill);
hz_20 = Fused_20/Fused_20(1);

[~, ~, Fused_50] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, CF50HzCa*C_Ca, SS, reserve_size, k_refill);
hz_50 = Fused_50/Fused_50(1);

%err_3KO = sqrt(sum((hz_1 - CFData3KO(1:20,1)).^2 + (hz_10 - CFData3KO(1:20,2)).^2 + (hz_20 - CFData3KO(1:20,3)).^2 + (hz_50(1:20) - CFData3KO(1:20,4)).^2) + 10*sum((hz_50(21:30) - CFData3KO(21:30,4)).^2));
err_3KO = sqrt(2*sum((hz_1(1:5) - CFData3KO(1:5,1)).^2 + (hz_10(1:5) - CFData3KO(1:5,2)).^2 + (hz_20(1:5) - CFData3KO(1:5,3)).^2 + (hz_50(1:5) - CFData3KO(1:5,4)).^2) + sum((hz_1(6:20) - CFData3KO(6:20,1)).^2 + (hz_10(6:20) - CFData3KO(6:20,2)).^2 + (hz_20(6:20) - CFData3KO(6:20,3)).^2 + (hz_50(6:20) - CFData3KO(6:20,4)).^2) + 10*sum((hz_50(21:30) - CFData3KO(21:30,4)).^2));

err = err_3KO;

%rec = [50 100 200 350 500 750 1000 2000 5000 10000];

% figure('Name','Syt3 KO Simulations','NumberTitle','off')
% subplot(3,2,1)
% plot(CFData3KO(1:20,1),'-k')
% title('1 Hz')
% xlabel('Pulse #')
% ylabel('Peak EPSC')
% set(gca,'xlim',[1 20])
% set(gca,'ylim',[0 1.2])
% hold on
% plot(hz_1,'ko','Markersize',5)
% %legend({'Data','Simulation'},'Position',[0.14, 0.75, .1, .1])
% 
% subplot(3,2,2)
% plot(CFData3KO(1:20,2),'-k')
% title('10 Hz')
% xlabel('Pulse #')
% ylabel('Peak EPSC')
% set(gca,'xlim',[1 20])
% set(gca,'ylim',[0 1])
% hold on
% plot(hz_10,'ko','Markersize',5)
% 
% subplot(3,2,3)
% plot(CFData3KO(1:20,3),'-k')
% title('20 Hz')
% xlabel('Pulse #')
% ylabel('Peak EPSC')
% set(gca,'xlim',[1 20])
% set(gca,'ylim',[0 1.2])
% hold on
% plot(hz_20,'ko','Markersize',5)
% 
% subplot(3,2,4)
% plot(CFData3KO(1:20,4),'-k')
% title('50 Hz')
% xlabel('Pulse #')
% ylabel('Peak EPSC')
% set(gca,'xlim',[1 20])
% set(gca,'ylim',[0 1.2])
% hold on
% plot(hz_50(1:20),'ko','Markersize',5)
% 
% subplot(3,2,5)
% semilogx(rec,CFData3KO(21:30,4),'-k')
% title('50 Hz Recovery')
% xlabel('t (ms)')
% ylabel('Peak EPSC')
% set(gca,'xlim',[50 10000])
% set(gca,'ylim',[0 1.2])
% hold on
% semilogx(rec,hz_50(21:30),'ko','Markersize',5)

disp(['x = [', num2str(x(1)), ',', num2str(x(2)), ',', num2str(x(3)), ',', num2str(x(4)), ',', num2str(x(5)), ',', num2str(x(6)), ']', ', err = ', num2str(err)])

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, Ca, SS, reserve_size, k_refill)

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
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, Ca(1,:), reserve_size, k_refill), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, Ca(i,:), reserve_size, k_refill), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,Ca,reserve_size,k_refill)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;

end

function dydt = dState(t,state,k_docking,k_undocking,Ca,reserve_size,k_refill)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;

end

end