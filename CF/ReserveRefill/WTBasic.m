CFDataWT = load('..\CFDataWT.mat').CFDataWT;

plot_data = 1;

%x = [0.70303,0.0033015,0.00010035,6.7483,0.0001,0]; %KO best fit 1/13/2021 err = 0.69845
%x = [0.87655,0.01415,0.0033386,2.7791,6.607e-05]; %50hz fit 3/1/21 cost = 0.417
x = [0.688,0.0030421,0.0032928,9.5651,6.607e-05]; %

p_release = x(1); 
k_docking = x(2); 
k_undocking = x(3); 
reserve_size = x(4);
k_refill = x(5);

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
[~,state] = ode113(@(t,state) dSS(t, state, k_docking, k_undocking, reserve_size, k_refill), [0 t_SS], state_0);

SS_WT = state(end,:);

[ts_1_WT, state_1_WT, Fused_1_WT] = stim_sim(stimulus_times_1, max_time_1, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill);
hz_1_WT = Fused_1_WT/Fused_1_WT(1);

[ts_10_WT, state_10_WT, Fused_10_WT] = stim_sim(stimulus_times_10, max_time_10, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill);
hz_10_WT = Fused_10_WT/Fused_10_WT(1);

[ts_20_WT, state_20_WT, Fused_20_WT] = stim_sim(stimulus_times_20, max_time_20, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill);
hz_20_WT = Fused_20_WT/Fused_20_WT(1);

[ts_50_WT, state_50_WT, Fused_50_WT] = stim_sim(stimulus_times_50, max_time_50, p_release, k_docking, k_undocking, SS_WT, reserve_size, k_refill);
hz_50_WT = Fused_50_WT/Fused_50_WT(1);

%err = sqrt(2*sum((hz_1_WT(1:5) - CFDataWT(1:5,1)).^2 + (hz_10_WT(1:5) - CFDataWT(1:5,2)).^2 + (hz_20_WT(1:5) - CFDataWT(1:5,3)).^2 + (hz_50_WT(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_1_WT(6:20) - CFDataWT(6:20,1)).^2 + (hz_10_WT(6:20) - CFDataWT(6:20,2)).^2 + (hz_20_WT(6:20) - CFDataWT(6:20,3)).^2 + (hz_50_WT(6:20) - CFDataWT(6:20,4)).^2) + 10*sum((hz_50_WT(21:30) - CFDataWT(21:30,4)).^2));
err = sqrt(2*sum((hz_50_WT(1:5) - CFDataWT(1:5,4)).^2) + sum((hz_50_WT(6:20) - CFDataWT(6:20,4)).^2));


disp(['Average error = ', num2str(err)])

% times = [0:0.1:ts_50_WT(end)]; %vector of length max_time denoting times from 1->max_time msec with step size 1
% peak = NaN(1,length(times));
% T_E = 2;
% alpha_1 = (times*exp(1)/T_E).*exp(-1*times/T_E); %reference alpha function
% alpha = zeros(1,length(times)); %functional alpha function
% 
% for i = 1:length(stimulus_times_50)
%     stimulus = stimulus_times_50(i);
%     peak((i-1)*200+1) = Fused_50_WT(i);
%     v = find(times == stimulus, 1, 'first'); %index of stimulus time
%     alpha(v:end) = alpha(v:end) + (Fused_50_WT(i)*alpha_1(1:end-v+1)); %calculate effect of each stimulus on the alpha function and sum them
% end
% figure('Name','EPSC Sim','NumberTitle','off')
% plot(times,alpha,'-r','LineWidth',2)
% hold on
% plot(times,peak,'ko','MarkerFaceColor','r','MarkerSize',20)
% set(gca,'Visible','off')
% title('Simulated 50Hz EPSC Train')
% xlabel('Time (ms)')
% ylabel('EPSC')
% set(gca,'xtick',[])
% set(gca,'ytick',[])

%export_fig('EPSC_Peak_Sim_50', '-dpng', '-transparent', '-nocrop', '-r300');
if plot_data == 1
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
end

function [ts, state, Fused] = stim_sim(stimulus_times, max_time, p_release, k_docking, k_undocking, SS, reserve_size, k_refill)

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
        [t,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end,:); out];

        ts = [ts(1:end,:); t];
    end
    
    if stimulus_times == linspace(0,20*19,20) %50hz
       
        Recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
        Fused_rec = zeros(length(Recovery),1);   
        
        for i = 1:length(Recovery)
            
            [~,out] = ode113(@(t,state) dState(t, state, k_docking, k_undocking, reserve_size, k_refill), stimulus_times(end)+[0 Recovery(i)], post_stim);
            pre_stim = out(end,:);
            Fused_rec(i) = pre_stim(2)*p_release;
       
        end
        Fused = [Fused; Fused_rec];
        
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
end

function dydt = dSS(~,state,k_docking,k_undocking,reserve_size,k_refill)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;

end

function dydt = dState(~,state,k_docking,k_undocking,reserve_size,k_refill)

dydt(1,1) = -state(1)*(state(3)/reserve_size)*k_docking + state(2)*k_undocking;
dydt(2,1) = -dydt(1,1);
dydt(3,1) = dydt(1,1) + (reserve_size-state(3))*k_refill;
end