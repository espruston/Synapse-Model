global p_mature
global p_immature
global k_docking
global k_undocking
global k_maturation
global k_dematuration
global delta_t
global SS
% global K_D_3
% global K_D_7
global k_on_3
global k_off_3
global C_3
global C_7
global Ca_rest
global Ca_spike
global Ca_residual
global T_Ca_decay
global sigma
global mu


%stimulus_times = [0 10];
stimulus_times = linspace(0,1000*19,20); %20 stims @ 1hz

% K_D_3 = 3e5; %M^-1ms^-1 Hui
% K_D_7 = 7.333e3; %M^-1ms^-1 Brandt/Knight 

k_on_3 = 3e5; %M^-1ms^-1 Hui
%k_on_3 = 3e4; %M^-1ms^-1
%k_on_3 = 0; %syt3 KO
k_off_3 = 0.05; %ms^-1  Hui

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

% hold values, PPR R^2 = .9724
% k_docking = .03;
% k_undocking = .0001;
% k_maturation = .000965;
% k_dematuration = .0001;

k_docking = .003;
k_undocking = 0;
k_maturation = .003;
k_dematuration = 0;

CDR = 0;
Facil = 0;

C_3 = 10;
C_7 = 0;

p_immature = .3;
p_mature = .7;

Ca_rest = 5e-8; %M
Ca_spike = 2e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms

delta_t = 1e-2; %ms

t_SS = 10000; %ms
ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [1; 0; 0; 0]; %start all vesicles in immature undocked state 

[t0,state] = ode15s(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,C_3,Ca_rest), [0 t_SS], state_0);

SS = state(end,:);

[Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50] = test4();

plot_four(Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50);

function [ts, state, Fused_im, Fused_m, Ca_sim] = stim_sim(stimulus_times, max_time)

    global p_immature
    global p_mature
    global k_docking
    global k_undocking
    global k_maturation
    global k_dematuration
    global k_on_3
    global k_off_3
    global C_3
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
        post_stim = pre_stim + [pre_stim(2)*p_immature+pre_stim(3)*p_mature, -pre_stim(2)*p_immature, -pre_stim(3)*p_mature, 0];
        Fused_im(i) = pre_stim(2)*p_immature;
        Fused_m(i) = pre_stim(3)*p_mature;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,C_3,Ca_sim), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end

    state = [SS; state];
    ts = [delta_t; ts];
end


function Ca_sim = create_Ca_signal(stimulus_times, max_time)
    
    global Ca_rest
    global Ca_spike
    global Ca_residual
    global T_Ca_decay
    global delta_t
    global sigma
    global mu
    
    
    ts = linspace(0,max_time,max_time/delta_t + 1);
    Ca_sim = zeros(1,length(ts));
    Ca_sim = Ca_sim + Ca_rest;

    for t = 1:length(stimulus_times) %simulate calcium influx

        spike_start_index = round(stimulus_times(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times(t)+mu)/delta_t) + 1;

        Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one

        Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

    end

end


function [Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50] = test4()
    
    stimulus_times_1 = linspace(0,1000*19,20); %20 stims 1hz
    stimulus_times_10 = linspace(0,100*19,20); %20 stims 10hz
    stimulus_times_20 = linspace(0,50*19,20); %20 stims 20hz
    stimulus_times_50 = linspace(0,20*19,20); 

    max_time_1 = stimulus_times_1(end) + stimulus_times_1(2)*3;
    [~, ~, Fused_im_1, Fused_m_1,~] = stim_sim(stimulus_times_1, max_time_1);
    Fused_1 = Fused_im_1 + Fused_m_1;
    hz_1 = Fused_1/Fused_1(1);
    
    max_time_10 = stimulus_times_10(end) + stimulus_times_10(2)*3;
    [~, ~, Fused_im_10, Fused_m_10,~] = stim_sim(stimulus_times_10, max_time_10);
    Fused_10 = Fused_im_10 + Fused_m_10;
    hz_10 = Fused_10/Fused_10(1);
    
    max_time_20 = stimulus_times_20(end) + stimulus_times_20(2)*3;
    [~, ~, Fused_im_20, Fused_m_20,~] = stim_sim(stimulus_times_20, max_time_20);
    Fused_20 = Fused_im_20 + Fused_m_20;
    hz_20 = Fused_20/Fused_20(1);

    max_time_50 = stimulus_times_50(end) + stimulus_times_50(2)*3;
    [~, ~, Fused_im_50, Fused_m_50, ~] = stim_sim(stimulus_times_50, max_time_50);
    Fused_50 = Fused_im_50 + Fused_m_50;
    hz_50 = Fused_50/Fused_50(1);

end

function plot_four(Fused_im_1, Fused_m_1, hz_1, Fused_im_10, Fused_m_10, hz_10, Fused_im_20, Fused_m_20, hz_20, Fused_im_50, Fused_m_50, hz_50)
    
    global C_3
    
    hz_1_data = [1 0.904585375 0.906090801 0.885396962 0.888286409 0.882153581 0.868080327 0.870438248 0.858248517 0.854839359 0.840352061 0.845592803 0.847296096 0.844539428 0.829714224 0.830002447 0.812125568 0.815760427 0.81009733 0.799518304];
    hz_10_data = [1 0.616949107 0.628326956 0.604464053 0.565850677 0.537103779 0.525309915 0.503692205 0.474574723 0.470769015 0.461639789 0.45151103 0.450073867 0.451884135 0.44312584 0.429802182 0.433756546 0.431020365 0.423930097 0.416302277];
    hz_20_data = [1 0.489671267 0.444366914 0.457655745 0.430402833 0.416929064 0.38834963 0.351531591 0.348392421 0.329990273 0.31629381 0.305756525 0.3022066 0.288966797 0.283932027 0.266221279 0.276873203 0.265228816 0.262831634 0.258863279];
    hz_50_data = [1 0.344586239 0.28753639 0.241585061 0.22098407 0.212199666 0.198432318 0.177849609 0.158522769 0.144448699 0.131731778 0.119584277 0.111943939 0.106093112 0.098933385 0.095593041 0.08975884 0.087444448 0.082325078 0.083824248];
    labels = ["1 Hz data", "10 Hz data", "20 Hz data", "50 Hz data"];
    
    if C_3 == 0 %KO
        hz_1_data = [1 0.883284703 0.858855554 0.845649686 0.837780118 0.842965361 0.828148214 0.806240914 0.799098827 0.792257423 0.789414849 0.785705274 0.763849243 0.724204838 0.750418701 0.740545466 0.743069941 0.74399981 0.736566288 0.736820291];
        hz_10_data = [1 0.452012818 0.44831423 0.433521032 0.376651907 0.361794245 0.322100309 0.326965448 0.309113437 0.303756475 0.293420616 0.298110805 0.286918883 0.27799293 0.292249154 0.292399085 0.284254389 0.286765184 0.277340747 0.278088675];
        hz_20_data = [1 0.348773535 0.295773884 0.279677667 0.268386245 0.235162278 0.21836187 0.203592063 0.18263844 0.175787404 0.156893656 0.154776094 0.149120775 0.1555285 0.142093573 0.138859878 0.132577603 0.136777797 0.132698411 0.126999323];
        hz_50_data = [1 0.240759282 0.216549504 0.161161579 0.129963837 0.117072816 0.108392052 0.097730189 0.090188357 0.084249752 0.081031452 0.077733613 0.072533983 0.075013212 0.071678789 0.069395147 0.066660259 0.065201118 0.066126873 0.066443286];
        labels = ["1 Hz KO data", "10 Hz KO data", "20 Hz KO data", "50 Hz KO data"];
        
    end
    
    Fused_1_norm = Fused_im_1(1) + Fused_m_1(1);
    Fused_10_norm = Fused_im_10(1) + Fused_m_10(1);
    Fused_20_norm = Fused_im_20(1) + Fused_m_20(1);
    Fused_50_norm = Fused_im_50(1) + Fused_m_50(1);
    
    figure 
    subplot(2,2,1)
    plot(hz_1_data,'-k')
    title('1 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 20])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_1,'ko','Markersize',5)
    plot(Fused_im_1/Fused_1_norm,'rv')
    plot(Fused_m_1/Fused_1_norm,'g^')
    %legend(labels(1),'1 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(2,2,2)
    plot(hz_10_data,'-k')
    title('10 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 20])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_10,'ko','Markersize',5)
    plot(Fused_im_10/Fused_10_norm,'rv')
    plot(Fused_m_10/Fused_10_norm,'g^')
    %legend(labels(2),'10 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(2,2,3)
    plot(hz_20_data,'-k')
    title('20 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 20])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_20,'ko','Markersize',5)
    plot(Fused_im_20/Fused_20_norm,'rv')
    plot(Fused_m_20/Fused_20_norm,'g^')
    %legend(labels(3),'20 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
    
    subplot(2,2,4)
    plot(hz_50_data,'-k')
    title('50 Hz')
    xlabel('Pulse #')
    ylabel('Peak EPSC')
    set(gca,'xlim',[1 20])
    set(gca,'ylim',[0 1])
    hold on
    plot(hz_50,'ko','Markersize',5)
    plot(Fused_im_50/Fused_50_norm,'rv')
    plot(Fused_m_50/Fused_50_norm,'g^')
    %legend(labels(4),'50 Hz Simulation','Low P Pool Fusion','High P Pool Fusion')
end

function plot_one(stimulus_times)
    
    global C_3
    
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

    hz_1_data = [1 0.904585375 0.906090801 0.885396962 0.888286409 0.882153581 0.868080327 0.870438248 0.858248517 0.854839359 0.840352061 0.845592803 0.847296096 0.844539428 0.829714224 0.830002447 0.812125568 0.815760427 0.81009733 0.799518304];
    hz_10_data = [1 0.616949107 0.628326956 0.604464053 0.565850677 0.537103779 0.525309915 0.503692205 0.474574723 0.470769015 0.461639789 0.45151103 0.450073867 0.451884135 0.44312584 0.429802182 0.433756546 0.431020365 0.423930097 0.416302277];
    hz_20_data = [1 0.489671267 0.444366914 0.457655745 0.430402833 0.416929064 0.38834963 0.351531591 0.348392421 0.329990273 0.31629381 0.305756525 0.3022066 0.288966797 0.283932027 0.266221279 0.276873203 0.265228816 0.262831634 0.258863279];
    hz_50_data = [1 0.344586239 0.28753639 0.241585061 0.22098407 0.212199666 0.198432318 0.177849609 0.158522769 0.144448699 0.131731778 0.119584277 0.111943939 0.106093112 0.098933385 0.095593041 0.08975884 0.087444448 0.082325078 0.083824248];

    if C_3 == 0 %KO
        hz_1_data = [1 0.883284703 0.858855554 0.845649686 0.837780118 0.842965361 0.828148214 0.806240914 0.799098827 0.792257423 0.789414849 0.785705274 0.763849243 0.724204838 0.750418701 0.740545466 0.743069941 0.74399981 0.736566288 0.736820291];
        hz_10_data = [1 0.452012818 0.44831423 0.433521032 0.376651907 0.361794245 0.322100309 0.326965448 0.309113437 0.303756475 0.293420616 0.298110805 0.286918883 0.27799293 0.292249154 0.292399085 0.284254389 0.286765184 0.277340747 0.278088675];
        hz_20_data = [1 0.348773535 0.295773884 0.279677667 0.268386245 0.235162278 0.21836187 0.203592063 0.18263844 0.175787404 0.156893656 0.154776094 0.149120775 0.1555285 0.142093573 0.138859878 0.132577603 0.136777797 0.132698411 0.126999323];
        hz_50_data = [1 0.240759282 0.216549504 0.161161579 0.129963837 0.117072816 0.108392052 0.097730189 0.090188357 0.084249752 0.081031452 0.077733613 0.072533983 0.075013212 0.071678789 0.069395147 0.066660259 0.065201118 0.066126873 0.066443286];

    end

    plot(hz_1_data,'r')
    plot(hz_10_data,'m')
    plot(hz_20_data,'g')
    plot(hz_50_data,'b')

    legend('Simulation','1 Hz data','10 Hz data','20 Hz data','50 Hz data')
    hold off
    
    subplot(3,2,2)
    semilogy(linspace(1,max_time,length(Ca_sim)),Ca_sim)
    ylabel('Simulated Calcium Concentration')
    xlabel('Time (ms)')
    set(gca,'xlim',[ts(1) max_time])

    subplot(3,2,3)
    hold on

    plot(ts, state(:,3),'-.b')
    plot(ts, state(:,4),'-b')
    plot(ts, state(:,3)+state(:,4),'-r')
    xlabel('Time (ms)')
    ylabel('Occupied Sites (norm)')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    hold off
    legend('Occupied Low P Sites','Occupied High P Sites','Total Occupied Sites');
    
    subplot(3,2,4)
    hold on

    plot(ts, state(:,1),'-.m')
    plot(ts, state(:,2),'-m')
    plot(ts, state(:,1)+state(:,2),'-r')
    xlabel('Time (ms)')
    ylabel('Empty Sites (norm)')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    hold off
    legend('Empty Low P Sites','Empty High P Sites','Total Empty Sites');

    subplot(3,2,5)
    plot(ts, state(:,5))
    xlabel('Time (ms)')
    ylabel('Syt3')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])

    subplot(3,2,6)
    plot(ts, state(:,6))
    xlabel('Time (ms)')
    ylabel('Syt7')
    set(gca,'xlim',[ts(1) max_time])
    set(gca,'ylim',[0 1])
    
end

function dydt = dSS(~,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,C_3,Ca_rest)
    
    k_maturation = k_maturation*(1+state(4)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    dydt(4,1) = (1-state(4))*k_on_3*Ca_rest - state(4)*k_off_3;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,C_3,Ca_sim)
    
    Ca = Ca_sim(round(t/.01)+1);
    k_maturation = k_maturation*(1+state(4)*C_3);
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    dydt(4,1) = (1-state(4))*k_on_3*Ca - state(4)*k_off_3;

end