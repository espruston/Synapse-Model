global p_mature
global p_immature
global k_docking
global k_undocking
global k_maturation
global k_dematuration
global delta_t
global SS
global k_on_3
global k_off_3
global k_on_7
global k_off_7
global C_3
global C_7
global Ca_rest
%global Ca_spike
global Ca_residual
global T_Ca_decay
%global sigma
global mu

type = "WT";

if type == "WT"
    data = matfile('WT_data.mat').WT_data;
elseif type == "syt3KO"
    data = matfile('syt3KO_data.mat').syt3KO_data;
elseif type == "syt7KO"
    data = matfile('syt7KO_data.mat').syt7KO_data;   
elseif type == "DKO"
    data = matfile('DKO_data.mat').DKO_data;
else
    disp("Type should be one of the following: WT, syt3KO, syt7KO, DKO \n Using WT")
    data = matfile('WT_data.mat').WT_data;
end

hz_1_data = data(:,1);
hz_10_data = data(:,2);
hz_20_data = data(:,3);
hz_50_data = data(:,4);
hz_100_data = data(:,5);
hz_200_data = data(:,6);
hz_200_recovery = data(1:11,7);

%stimulus_times = [0 10];
stimulus_times = linspace(0,5*99,100); %100 stims @ 200hz

k_on_3 = 3e5; %M^-1ms^-1 Hui
k_off_3 = 0.05; %ms^-1  Hui
k_on_7 = 7.333e3; %Knight
k_off_7 = 1.1e-2;

FWHM = .34; %Local calcium full width half maximum ms
%sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

% hold values, PPR R^2 = .9724
% k_docking = .03;
% k_undocking = .0001;
% k_maturation = .000965;
% k_dematuration = .0001;

k_docking = 0.0033;
k_undocking = 0.00039;
k_maturation = .00019;
k_dematuration = .000045;

p_immature = .15;
p_mature = .53;

C_3 = 1.3;
C_7 = 1;

Ca_rest = 50e-9; %M
%Ca_spike = 2e-5; %M
Ca_residual = 500e-9; %M
T_Ca_decay = 50; %ms

delta_t = 1e-2; %ms

t_SS = 10000; %ms
ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [1; 0; 0; 0; 0]; %[empty; immature; mature; syt3; syt7]

[t0,state] = ode15s(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_rest), [0 t_SS], state_0);

SS = state(end,:);

[Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50, Fused_im_100, Fused_m_100, hz_100, Fused_im_200, Fused_m_200, hz_200] = test6();

plot_seven(data,Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50, Fused_im_100, Fused_m_100, hz_100, Fused_im_200, Fused_m_200, hz_200);

%plot_one(data,stimulus_times)
    
function [ts, state, Fused_im, Fused_m, Ca_sim] = stim_sim(stimulus_times, max_time)

    global p_immature
    global p_mature
    global k_docking
    global k_undocking
    global k_maturation
    global k_dematuration
    global k_on_3
    global k_off_3
    global k_on_7
    global k_off_7
    global C_3
    global C_7
    global delta_t
    global SS
    
    state = SS;
    Ca_sim = create_Ca_signal(stimulus_times, max_time);
    %[syt3, syt7] = syt_sim(Ca_sim);
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused_im = zeros(length(stimulus_times),1);
    Fused_m = zeros(length(stimulus_times),1);

    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_immature*(1+C_7*pre_stim(5))+pre_stim(3)*p_mature*(1+C_7*pre_stim(5)), -pre_stim(2)*p_immature*(1+C_7*pre_stim(5)), -pre_stim(3)*p_mature*(1+C_7*pre_stim(5)), 0, 0];
        Fused_im(i) = pre_stim(2)*p_immature*(1+C_7*pre_stim(5));
        Fused_m(i) = pre_stim(3)*p_mature*(1+C_7*pre_stim(5));
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_sim), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end

    state = [SS; state];
    ts = [delta_t; ts];
end

%only simulate res. Ca
function Ca_sim = create_Ca_signal(stimulus_times, max_time)
    
    global Ca_rest
    %global Ca_spike
    global Ca_residual
    global T_Ca_decay
    global delta_t
    %global sigma
    global mu
    
    
    ts = linspace(0,max_time,max_time/delta_t + 1);
    Ca_sim = zeros(1,length(ts));
    Ca_sim = Ca_sim + Ca_rest;

    for t = 1:length(stimulus_times) %simulate calcium influx

        %spike_start_index = round(stimulus_times(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times(t)+mu)/delta_t) + 1;

        %Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one

        Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

    end

end


function [Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50, Fused_im_100, Fused_m_100, hz_100, Fused_im_200, Fused_m_200, hz_200] = test6()
    
    stimulus_times_1 = linspace(0,1000*99,100); %100 stims 1hz
    stimulus_times_10 = linspace(0,100*99,100); %100 stims 10hz
    stimulus_times_20 = linspace(0,50*99,100); 
    stimulus_times_50 = linspace(0,20*99,100);
    stimulus_times_100 = linspace(0,10*99,100);
    stimulus_times_200 = linspace(0,5*99,100);
    stimulus_times_200 = [stimulus_times_200, stimulus_times_200(end)+10, stimulus_times_200(end)+20, stimulus_times_200(end)+50, stimulus_times_200(end)+100, stimulus_times_200(end)+200, stimulus_times_200(end)+500, stimulus_times_200(end)+1000, stimulus_times_200(end)+2000, stimulus_times_200(end)+5000, stimulus_times_200(end)+10000, stimulus_times_200(end)+20000]; %recovery


    max_time_1 = stimulus_times_1(end) + stimulus_times_1(2)*30;
    [~, ~, Fused_im_1, Fused_m_1,~] = stim_sim(stimulus_times_1, max_time_1);
    Fused_1 = Fused_im_1 + Fused_m_1;
    hz_1 = Fused_1/Fused_1(1);
    
    max_time_10 = stimulus_times_10(end) + stimulus_times_10(2)*30;
    [~, ~, Fused_im_10, Fused_m_10,~] = stim_sim(stimulus_times_10, max_time_10);
    Fused_10 = Fused_im_10 + Fused_m_10;
    hz_10 = Fused_10/Fused_10(1);
    
    max_time_20 = stimulus_times_20(end) + stimulus_times_20(2)*30;
    [~, ~, Fused_im_20, Fused_m_20,~] = stim_sim(stimulus_times_20, max_time_20);
    Fused_20 = Fused_im_20 + Fused_m_20;
    hz_20 = Fused_20/Fused_20(1);

    max_time_50 = stimulus_times_50(end) + stimulus_times_50(2)*30;
    [~, ~, Fused_im_50, Fused_m_50, ~] = stim_sim(stimulus_times_50, max_time_50);
    Fused_50 = Fused_im_50 + Fused_m_50;
    hz_50 = Fused_50/Fused_50(1);
    
    max_time_100 = stimulus_times_100(end) + stimulus_times_100(2)*100;
    [~, ~, Fused_im_100, Fused_m_100, ~] = stim_sim(stimulus_times_100, max_time_100);
    Fused_100 = Fused_im_100 + Fused_m_100;
    hz_100 = Fused_100/Fused_100(1);
    
    max_time_200 = stimulus_times_200(end) + stimulus_times_200(2)*30;
    [~, ~, Fused_im_200, Fused_m_200, ~] = stim_sim(stimulus_times_200, max_time_200);
    Fused_200 = Fused_im_200 + Fused_m_200;
    hz_200 = Fused_200/Fused_200(1);

end


function plot_seven(data,Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50, Fused_im_100, Fused_m_100, hz_100, Fused_im_200, Fused_m_200, hz_200)
    
    hz_1_data = data(:,1);
    hz_10_data = data(:,2);
    hz_20_data = data(:,3);
    hz_50_data = data(:,4);
    hz_100_data = data(:,5);
    hz_200_data = data(:,6);
    hz_200_recovery = data(1:11,7);

    labels = ["1 Hz data", "10 Hz data", "20 Hz data", "50 Hz data", "100 Hz data", "200 Hz data"];
    
    Fused_1_norm = Fused_im_1(1) + Fused_m_1(1);
    Fused_10_norm = Fused_im_10(1) + Fused_m_10(1);
    Fused_20_norm = Fused_im_20(1) + Fused_m_20(1);
    Fused_50_norm = Fused_im_50(1) + Fused_m_50(1);
    Fused_100_norm = Fused_im_100(1) + Fused_m_100(1);
    Fused_200_norm = Fused_im_200(1) + Fused_m_200(1);
    
    figure 
    subplot(4,2,1)
    plot(hz_1_data,'-k')
    title('1 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_1,'ko','Markersize',5)
    plot(Fused_im_1/Fused_1_norm,'rv')
    plot(Fused_m_1/Fused_1_norm,'g^')
    legend(labels(1),'1 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(4,2,2)
    plot(hz_10_data,'-k')
    title('10 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_10,'ko','Markersize',5)
    plot(Fused_im_10/Fused_10_norm,'rv')
    plot(Fused_m_10/Fused_10_norm,'g^')
    %legend(labels(2),'10 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(4,2,3)
    plot(hz_20_data,'-k')
    title('20 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_20,'ko','Markersize',5)
    plot(Fused_im_20/Fused_20_norm,'rv')
    plot(Fused_m_20/Fused_20_norm,'g^')
    %legend(labels(3),'20 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(4,2,4)
    plot(hz_50_data,'-k')
    title('50 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_50,'ko','Markersize',5)
    plot(Fused_im_50/Fused_50_norm,'rv')
    plot(Fused_m_50/Fused_50_norm,'g^')
    %legend(labels(4),'50 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(4,2,5)
    plot(hz_100_data,'-k')
    title('100 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_100,'ko','Markersize',5)
    plot(Fused_im_100/Fused_100_norm,'rv')
    plot(Fused_m_100/Fused_100_norm,'g^')
    
    subplot(4,2,6)
    plot(hz_200_data(1:100),'-k')
    title('200 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 100])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_200,'ko','Markersize',5)
    plot(Fused_im_200/Fused_200_norm,'rv')
    plot(Fused_m_200/Fused_200_norm,'g^')
    
    subplot(4,2,7)
    semilogx([10,20,50,100,200,500,1000,2000,5000,10000,20000],hz_200_recovery,'-k')
    title('200 Hz Recovery')
    xlabel('Time after Last Pulse (ms)')
    ylabel('EPSC (norm)')
    set(gca,'ylim',[0 1.2])
    hold on
    semilogx([10,20,50,100,200,500,1000,2000,5000,10000,20000],hz_200(101:111),'ko','Markersize',5)
    semilogx([10,20,50,100,200,500,1000,2000,5000,10000,20000],Fused_im_200(101:111)/Fused_200_norm,'rv')
    semilogx([10,20,50,100,200,500,1000,2000,5000,10000,20000],Fused_m_200(101:111)/Fused_200_norm,'g^')
end

function plot_one(data, stimulus_times)

    max_time = stimulus_times(end) + stimulus_times(2)*3;
    [ts, state, Fused_im, Fused_m, Ca_sim] = stim_sim(stimulus_times, max_time);
    Fused = Fused_im + Fused_m;
    Fused = Fused/Fused(1);
    
    figure
    subplot(3,2,1)
    plot(Fused,'ko','MarkerSize', 5)
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 length(Fused)])
    set(gca,'ylim',[0 1])
    hold on
    
    hz_1_data = data(:,1);
    hz_10_data = data(:,2);
    hz_20_data = data(:,3);
    hz_50_data = data(:,4);
    hz_100_data = data(:,5);
    hz_200_data = data(:,6);

    labels = ["1 Hz data", "10 Hz data", "20 Hz data", "50 Hz data", "100 Hz data", "200 Hz data"];
    

    plot(hz_1_data,'b')
    plot(hz_10_data,'Color',[17 17 17]/255)
    plot(hz_20_data,'g')
    plot(hz_50_data,'m')
    plot(hz_100_data,'r')
    plot(hz_200_data,'y')

    legend('Simulation','1 Hz data','10 Hz data','20 Hz data','50 Hz data','100 Hz data','200 Hz data')
    hold off
    
    subplot(3,2,2)
    semilogy(linspace(1,max_time,length(Ca_sim)),Ca_sim)
    ylabel('Simulated Calcium Concentration')
    xlabel('Time (ms)')
    set(gca,'xlim',[ts(1) max_time])

    subplot(3,2,3)
    hold on

    plot(ts, state(:,2),'-.b')
    plot(ts, state(:,3),'-b')
    plot(ts, state(:,2)+state(:,3),'-r')
    xlabel('Time (ms)')
    ylabel('Occupied Sites (norm)')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    hold off
    legend('Occupied Low P Sites','Occupied High P Sites','Total Occupied Sites');
    
    subplot(3,2,4)
    hold on

    plot(ts, state(:,1),'-r')
    xlabel('Time (ms)')
    ylabel('Empty Sites (norm)')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    hold off

    subplot(3,2,5)
    plot(ts, state(:,4))
    xlabel('Time (ms)')
    ylabel('Syt3')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    subplot(3,2,6)
    plot(ts, state(:,5))
    xlabel('Time (ms)')
    ylabel('Syt7')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])
    
end


function dydt = dSS(~,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_rest)
    
    k_docking = k_docking*(1+state(4)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    dydt(4,1) = (1-state(4))*k_on_3*Ca_rest - state(4)*k_off_3;
    dydt(5,1) = (1-state(5))*k_on_7*Ca_rest - state(5)*k_off_7;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_sim)
    
    Ca = Ca_sim(round(t/.01)+1);
    k_docking = k_docking*(1+state(4)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    dydt(4,1) = (1-state(4))*k_on_3*Ca - state(4)*k_off_3;
    dydt(5,1) = (1-state(5))*k_on_7*Ca - state(5)*k_off_7;
    
end
