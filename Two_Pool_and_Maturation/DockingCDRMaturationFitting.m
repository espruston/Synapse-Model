%stimulus_times = [0 10];
%stimulus_times = linspace(0,20*19,20); %20 stims @50hz

global K_D_3
global K_D_7

K_D_3 = 7e-6; %1/2Ca M syt3
K_D_7 = 5.4e-6; %1/2Ca M syt7

% hold values, PPR R^2 = .93419, SS R^2 = .97901 train recovery R^2 = .97743, combined R^2 = .97823
% k_docking = .01;
% k_undocking = 0;
% k_maturation = 0.003;
% k_dematuration = 0;
% 
% C_3 = 50; %concentration of syt3
% C_7 = 0; %concentration of syt7
% 
% CDR = .001; %effect of calcium dependent recovery from depression
% Facil = 0; %effect of facilitation
% 
% p_immature = .1;
% p_mature = .8;

global k_docking
global k_undocking
global k_maturation
global k_dematuration

global C_3
global C_7

global CDR
global Facil

global p_immature
global p_mature

k_docking = .01;
k_undocking = 0;
k_maturation = 0.003;
k_dematuration = 0;

C_3 = 50; %concentration of syt3
C_7 = 0; %concentration of syt7

CDR = .001; %effect of calcium dependent recovery from depression
Facil = 0; %effect of facilitation

p_immature = .1;
p_mature = .8;


global Ca_rest
global Ca_spike
global Ca_residual
global T_Ca_decay

global sigma
global mu

Ca_rest = 5e-8; %M
Ca_spike = 2e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

global delta_t
delta_t = 1e-2; %ms

t_SS = 100000; %ms
ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [0; 1; 0]; %start all vesicles in immature docked state 

%[~,state] = ode45(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration), [0 t_SS], state_0);
[~,state] = ode45(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration,(Ca_rest^2/(K_D_3^2 + Ca_rest^2))*C_3, CDR), [0 t_SS], state_0);

global SS
SS = state(end,:);

[ISI_PPR, PPR, R_PPR, PPR_data, PPR_data_SEMs] = PPR_fit();

[freq, trainSS, R_trainSS, trainSS_data] = trainSS_fit();

[ISI_train_recovery, train_recovery, R_train_recovery, train_recovery_data, train_recovery_data_SEMs] = trainRecovery_fit();

subplot(3,1,1)
semilogx(ISI_PPR,PPR,'-bo')
title('PPR simulation and data')
xlabel('Time(ms)')
ylabel('PPR')
hold on
errorbar(ISI_PPR,PPR_data,PPR_data_SEMs,'-r')
legend('Simulation','Data')
hold off

subplot(3,1,2)
semilogx(ISI_train_recovery,train_recovery,'-bo')
title('Train recovery simulation and data')
xlabel('Time(ms)')
ylabel('Recovery')
hold on
errorbar(ISI_train_recovery,train_recovery_data,train_recovery_data_SEMs,'-r')
legend('Simulation','Data')
hold off

subplot(3,1,3)
plot(freq,trainSS,'-bo')
title('Train steady state simulation and data')
xlabel('Frequency (Hz)')
ylabel('Steady state EPSC')
hold on
plot(freq,trainSS_data,'-rx')
legend('Simulation','Data')
hold off

R_squared = [R_PPR.^2 R_train_recovery.^2 R_trainSS.^2];
R_squared = transpose(R_squared([2 6 10]));
R_squared = [R_squared; mean(R_squared)];
Test = {'PPR'; 'Train Recovery'; 'Train Steady State'; 'Combined'};
disp(table(Test, R_squared))

%figure
% subplot(2,2,1)
% plot(ts, [0; 0; EPSC_sim])
% xlabel('Time (ms)')
% ylabel('EPSC (norm)')
% set(gca,'xlim',[ts(1) max_time])
% 
% subplot(2,2,2)
% plot(ts, state(:,1))
% xlabel('Time (ms)')
% ylabel('Empty sites')
% set(gca,'xlim',[ts(1) max_time])
% 
% subplot(2,2,3)
% plot(ts, state(:,2))
% xlabel('Time (ms)')
% ylabel('Docked sites')
% set(gca,'xlim',[ts(1) max_time])
% 
% subplot(2,2,4)
% plot(ts, state(:,3))
% xlabel('Time (ms)')
% ylabel('Mature sites')
% set(gca,'xlim',[ts(1) max_time])


function [ts, state, Fused, Ca_sim, syt3, syt7] = stim_sim(stimulus_times, max_time)

    global p_immature
    global p_mature
    global k_docking
    global k_undocking
    global k_maturation
    global k_dematuration
    global SS
    global CDR
    
    state = SS;
    Ca_sim = create_Ca_signal(stimulus_times, max_time);
    [syt3, syt7] = syt_sim(Ca_sim);
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end);

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_immature + pre_stim(3)*p_mature, -pre_stim(2)*p_immature, -pre_stim(3)*p_mature];
        Fused(i) = pre_stim(2)*p_immature + pre_stim(3)*p_mature;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration, syt3, CDR), [0 stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1); t+ts(end)];
    end

end


function [ISI, PPR, R_PPR, PPR_data, PPR_data_SEMs] = PPR_fit()
    
    ISI = [10 20 35 50 75 100 200 350 500 750 1000 2000 3500 5000 7500 10000];
    PPR = zeros(1,length(ISI));

    for j = 1:length(ISI)
        stimulus_times = [0 ISI(j)];

        max_time = stimulus_times(end) + stimulus_times(2)*3;
        
        [~, ~, Fused, ~, ~] = stim_sim(stimulus_times, max_time);
        Fused = Fused/Fused(1);
        PPR(j) = Fused(2);

    end

    PPR_data_WT = [0.1981008 0.298926974 0.361852919 0.389062492 0.468679125 0.559630866 0.740901035 0.806576031 0.840212316 0.878647834 0.907643109 0.965122849 0.96302892	0.990415868	1.007491546	0.999443074];
    PPR_data_SEMs_WT = [0.025540668	0.038168232	0.044923775	0.046957946	0.045900835	0.043358642	0.024948038	0.027493994	0.023202978	0.018398629	0.01439941 0.00835564 0.017361754 0.01835476 0.01216081	0.007243546];
    
    %PPR_data_WT_Ca_block = [.103 .194 .272 .307 .346 .365 .417 .477 .503 .620 .615 .850 .936 .965 .980 1.04]; 
    %PPR_data_SEMs_WT_Ca_block = [0.023707394 0.031631686 0.036341527 0.024351258 0.020447046 0.022314029 0.02242531	0.024670909	0.024094207	0.020289501	0.019038797	0.012231847	0.012046944	0.006259344	0.010283128	0.00864487];

    %PPR_data_syt3KO_Ca_block = [0.151147669	0.184327088	0.238073419	0.278179684	0.349569587	0.371625114	0.4931594 0.601356914 0.617759941 0.676810496 0.743110597 0.876359759 0.944711679 0.973497529 0.976220011 1.017493162];
    %PPR_data_SEMs_syt3KO_Ca_block = [0.026678818 0.023828808 0.033945846 0.031585123 0.030657408 0.032225669 0.022723458 0.051859004 0.053610984 0.051481764 0.040548383 0.028425046 0.021286098 0.018492625 0.01164004	0.010175575];
    
    %PPR_data = PPR_data_syt3KO_Ca_block;
    %PPR_data_SEMs = PPR_data_SEMs_syt3KO_Ca_block;
    
    PPR_data = PPR_data_WT;
    PPR_data_SEMs = PPR_data_SEMs_WT;
    
    R_PPR = corrcoef(PPR, PPR_data);

end


function [freq, trainSS, R_trainSS, trainSS_data] = trainSS_fit()
    
    freq = [1 10 20 50];
    trainSS = zeros(1,length(freq));

    for j = 1:length(freq)
        stimulus_times = linspace(0, 1000/freq(j)*19, 20);

        max_time = stimulus_times(end) + stimulus_times(2)*10;
        
        [~, ~, Fused, ~, ~] = stim_sim(stimulus_times, max_time);

        Fused = Fused/Fused(1);
        trainSS(j) = mean(Fused(end-5:end));

    end

    trainSS_data_WT = [0.835151796 0.396651642 0.214972754 0.008022113];
    
    %trainSS_data_WT_Ca_block = [.722 .241 .104 .0138]; 

    %trainSS_data_syt3KO_Ca_block = [0.761628403 0.276455062 0.127877809 0.016273454];
    
    %trainSS_data = trainSS_data_syt3KO_Ca_block;
    
    trainSS_data = trainSS_data_WT;
    
    R_trainSS = corrcoef(trainSS, trainSS_data);

end


function [ISI, recovery, R_train_recovery, train_recovery_data, train_recovery_data_SEMs] = trainRecovery_fit()
    
    global p_immature
    global p_mature
    global k_docking
    global k_undocking
    global k_maturation
    global k_dematuration
    global CDR
    
    stimulus_times = linspace(0,20*19,20); %20 stims @50hz

    max_time = stimulus_times(end) + stimulus_times(2)*10;
    
    [ts, state, Fused, ~, syt3, ~] = stim_sim(stimulus_times, max_time);

    end_index = find(ts==stimulus_times(end));
    post_train = state(end_index,:);
    syt3 = syt3(end_index:end);
    ISI = [50 100 200 350 500 750 1000 2000 5000 10000];
    recovery = zeros(1,length(ISI));

    for j = 1:length(ISI)

        [~,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,syt3,CDR), [0 ISI(j)], post_train);
        recovery(j) = (out(end,2)*p_immature + out(end,3)*p_mature)/Fused(1);

    end

    train_recovery_data_WT = [0.102997317 0.241270486 0.528117985 0.665949583 0.807115869 0.883204702 0.869620812 0.978624258 1.009330864 1.016249054];
    train_recovery_data_SEMs_WT = [0.032003424 0.045835464 0.070158349 0.07488678 0.074286201 0.070354653 0.071985604 0.025662978 0.037203997 0.031382848];
    
    %train_recovery_data_WT_Ca_block = [.133 .301 .535 .706 .853 .918 .950 1.03 1.01 1.03];
    %train_recovery_data_SEMs_WT_Ca_block = [0.011220544	0.02984717 0.052266679 0.048872772 0.029092336 0.026307037 0.021643624 0.012835957 0.009772508 0.027547534];
    
    %train_recovery_data_syt3KO_Ca_block = [0.095687504 0.1867 0.483622448 0.832322146 1.02638655 1.15260018	1.049615847	0.890841363	1.010482505	1.002926784];
    %train_recovery_data_SEMs_syt3KO_Ca_block = [0.020921145	0.025491913	0.053612745	0.036642093	0.023708192	0.058597074	0.020695876	0.066537454	0.012680184	0.017007916];
    
    %train_recovery_data = train_recovery_data_syt3KO_Ca_block;
    %train_recovery_data_SEMs = train_recovery_data_SEMs_syt3KO_Ca_block;
    
    train_recovery_data = train_recovery_data_WT;
    train_recovery_data_SEMs = train_recovery_data_SEMs_WT;
    R_train_recovery = corrcoef(recovery, train_recovery_data);
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


function [syt3, syt7] = syt_sim(Ca_sim)
    
    global K_D_3
    global K_D_7
    global C_3
    global C_7
    
    Ca_squared = Ca_sim.^2;
    
    syt3 = (Ca_squared./(K_D_3^2 + Ca_squared)).*C_3;
    syt7 = (Ca_squared./(K_D_7^2 + Ca_squared)).*C_7;

end


function dydt = dSS(~,state,k_docking,k_undocking,k_maturation,k_dematuration,syt3,CDR)
    
    dydt(1) = - state(1)*k_docking + state(2)*k_undocking;
    dydt(2) = state(1)*k_docking - state(2)*k_undocking - state(2)*(k_maturation + syt3*CDR) + state(3)*k_dematuration;
    dydt(3) = state(2)*(k_maturation + syt3*CDR) - state(3)*k_dematuration;
    dydt = transpose(dydt);

end


function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,syt3,CDR)
    
    syt3 = syt3(round(t)+1);
    
    dydt(1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2) = state(1)*k_docking - state(2)*k_undocking - state(2)*(k_maturation + syt3*CDR) + state(3)*k_dematuration;
    dydt(3) = state(2)*(k_maturation + syt3*CDR) - state(3)*k_dematuration;
    dydt = transpose(dydt);
end