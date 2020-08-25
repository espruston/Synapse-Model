global Ca_rest
global Ca_spike
global Ca_residual
global T_Ca_decay
global delta_t
global sigma
global mu
global p_release
global k_docking
global k_undocking
global k_docking_transient
global k_undocking_transient
global SS
global K_D_3
global K_D_7
global C_3
global C_7

Ca_rest = 5e-8; %M
Ca_spike = 2e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms
delta_t = .01;

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)


p_release = .8;

k_docking = .0016;
k_undocking = .0017;

k_docking_transient = 0.01;
k_undocking_transient = 0.005;

K_D_3 = 7e-6; %1/2Ca M syt3
K_D_7 = 5.4e-6; %1/2Ca M syt7

C_3 = 0; %concentration of syt3
C_7 = 0; %concentration of syt7

% 
% k_docking_transient = 0;
% k_undocking_transient = 0;

t_SS = 10000;

state_0 = [1; 0]; %all sites start empty
[~,state] = ode45(@(t,state) dSS(t,state,k_docking,k_undocking,k_docking_transient,k_undocking_transient,0,0,Ca_rest), [0 t_SS], state_0);

SS = state(end,:);

ISI_PPR = [10 20 35 50 75 100 200 350 500 750 1000 2000 3500 5000 7500 10000];
PPR = zeros(1,length(ISI_PPR));

for j = 1:length(ISI_PPR)
    stimulus_times = [0 ISI_PPR(j)];

    max_time = stimulus_times(end) + stimulus_times(2)*3;

    [~, ~, Fused] = stim_sim(stimulus_times, max_time);
    Fused = Fused/Fused(1);
    PPR(j) = Fused(2);

end
% 
PPR_data_WT = [0.22 0.32 0.39 0.42 0.50 0.58 0.70 0.75 0.80 0.84 0.87 0.94 0.97 0.98 1.00 1.01]; %set 8/23/20
PPR_data_SEMs_WT = [0.02 0.03 0.03 0.03 0.02 0.02 0.03 0.03 0.02 0.02 0.02 0.01 0.01 0.02 0.01 0.01]; %set 8/23/20

% PPR_data_WT_Ca_block = [.103 .194 .272 .307 .346 .365 .417 .477 .503 .620 .615 .850 .936 .965 .980 1.04]; 
% PPR_data_SEMs_WT_Ca_block = [0.023707394 0.031631686 0.036341527 0.024351258 0.020447046 0.022314029 0.02242531	0.024670909	0.024094207	0.020289501	0.019038797	0.012231847	0.012046944	0.006259344	0.010283128	0.00864487];

% PPR_data_syt3KO_Ca_block = [0.151147669	0.184327088	0.238073419	0.278179684	0.349569587	0.371625114	0.4931594 0.601356914 0.617759941 0.676810496 0.743110597 0.876359759 0.944711679 0.973497529 0.976220011 1.017493162];
% PPR_data_SEMs_syt3KO_Ca_block = [0.026678818 0.023828808 0.033945846 0.031585123 0.030657408 0.032225669 0.022723458 0.051859004 0.053610984 0.051481764 0.040548383 0.028425046 0.021286098 0.018492625 0.01164004	0.010175575];
% 
% PPR_data = PPR_data_syt3KO_Ca_block;
% PPR_data_SEMs = PPR_data_SEMs_syt3KO_Ca_block;

PPR_data = PPR_data_WT;
PPR_data_SEMs = PPR_data_SEMs_WT;

% PPR_data = PPR_data_WT_Ca_block;
% PPR_data_SEMs = PPR_data_SEMs_WT_Ca_block;

R_PPR = corrcoef(PPR, PPR_data);


stimulus_times = linspace(0,20*19,20); %20 stims @50hz

max_time = stimulus_times(end) + stimulus_times(2)*10;

[ts, state, ~] = stim_sim(stimulus_times, max_time);

post_train = state(ts==stimulus_times(end),:);
Ca = create_Ca_signal(stimulus_times, stimulus_times(end)+10000);

ISI_train_recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
syt3 = zeros(10001,1);
syt7 = zeros(10001,1);
train_recovery = zeros(1,length(ISI_train_recovery));

for j = 1:length(ISI_train_recovery)

    [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_docking_transient,k_undocking_transient,syt3,syt7,Ca(stimulus_times(end)/delta_t:end)/Ca_spike), [0 ISI_train_recovery(j)], post_train);
    train_recovery(j) = out(end,2)*p_release/Fused(1);

end

train_recovery_data_WT = [0.214594108 0.384408264 0.638176057 0.817518904 0.887058714 0.872159256 0.964188913 0.963314927 0.982736016 0.981537505]; %set 8/23/20
train_recovery_data_SEMs_WT = [0.035160387 0.040212604 0.043887001 0.027611691 0.025484558 0.05457874 0.03331073 0.024722031 0.011843753 0.011107924]; %set 8/23/20

% train_recovery_data_WT_Ca_block = [.133 .301 .535 .706 .853 .918 .950 1.03 1.01 1.03];
% train_recovery_data_SEMs_WT_Ca_block = [0.011220544	0.02984717 0.052266679 0.048872772 0.029092336 0.026307037 0.021643624 0.012835957 0.009772508 0.027547534];

% train_recovery_data_syt3KO_Ca_block = [0.095687504 0.1867 0.483622448 0.832322146 1.02638655 1.15260018	1.049615847	0.890841363	1.010482505	1.002926784];
% train_recovery_data_SEMs_syt3KO_Ca_block = [0.020921145	0.025491913	0.053612745	0.036642093	0.023708192	0.058597074	0.020695876	0.066537454	0.012680184	0.017007916];
% 
% train_recovery_data = train_recovery_data_syt3KO_Ca_block;
% train_recovery_data_SEMs = train_recovery_data_SEMs_syt3KO_Ca_block;

train_recovery_data = train_recovery_data_WT;
train_recovery_data_SEMs = train_recovery_data_SEMs_WT;

% train_recovery_data = train_recovery_data_WT_Ca_block;
% train_recovery_data_SEMs = train_recovery_data_SEMs_WT_Ca_block;

R_train_recovery = corrcoef(train_recovery, train_recovery_data);


freq = [1 10 20 50];
trainSS = zeros(1,length(freq));

for j = 1:length(freq)
    
    stimulus_times = linspace(0, 1000/freq(j)*19, 20);

    max_time = stimulus_times(end) + stimulus_times(2)*10;

    [~, ~, Fused] = stim_sim(stimulus_times, max_time);

    Fused = Fused/Fused(1);
    trainSS(j) = mean(Fused(end-5:end));

end

trainSS_data_WT = [0.835151796 0.396651642 0.2587 0.008022113]; %20hz freq set 8/23/20
trainSS_data_SEMs_WT = [0.013 0.00933 0.0077 0.00357]; %20hz set 8/23/20

% trainSS_data_WT_Ca_block = [.722 .241 .104 .0138]; 

% trainSS_data_syt3KO_Ca_block = [0.761628403 0.276455062 0.127877809 0.016273454];
% 
% trainSS_data = trainSS_data_syt3KO_Ca_block;

trainSS_data = trainSS_data_WT;
trainSS_data_SEMs = trainSS_data_SEMs_WT;

% trainSS_data = trainSS_data_WT_Ca_block;

R_trainSS = corrcoef(trainSS, trainSS_data);

subplot(3,1,1)
semilogx(ISI_PPR,PPR,'-b.')
title('PPR simulation and data')
xlabel('Time(ms)')
ylabel('PPR')
hold on
errorbar(ISI_PPR,PPR_data,PPR_data_SEMs,'-r')
legend('Simulation','Data')
hold off

subplot(3,1,2)
semilogx(ISI_train_recovery,train_recovery,'-b.')
title('Train recovery simulation and data')
xlabel('Time(ms)')
ylabel('Recovery')
hold on
errorbar(ISI_train_recovery,train_recovery_data,train_recovery_data_SEMs,'-r')
legend('Simulation','Data')
hold off

subplot(3,1,3)
plot(freq,trainSS,'-b.')
title('Train steady state simulation and data')
xlabel('Frequency (Hz)')
ylabel('Steady state EPSC')
hold on
errorbar(freq,trainSS_data,trainSS_data_SEMs,'-r')
legend('Simulation','Data')
hold off


sq_sum_PPR = sum(((PPR-PPR_data_WT)./PPR_data_SEMs_WT./16).^2);
sq_sum_train_recovery = sum(((train_recovery-train_recovery_data_WT)./train_recovery_data_SEMs_WT./10).^2);
sq_sum_trainSS = sum(((trainSS-trainSS_data_WT)./trainSS_data_SEMs_WT./4).^2);

sq_sum_total =  sq_sum_trainSS + sq_sum_PPR + sq_sum_train_recovery;


disp(sq_sum_trainSS/sq_sum_total)
disp(sq_sum_train_recovery/sq_sum_total)
disp(sq_sum_PPR/sq_sum_total)


function [ts, state, Fused] = stim_sim(stimulus_times, max_time)

    global p_release
    global k_docking
    global k_undocking
    global k_docking_transient
    global k_undocking_transient
    global delta_t
    global SS
    global Ca_spike
    
    state = SS;
    Ca_sim = create_Ca_signal(stimulus_times, max_time);
    [syt3, syt7] = syt_sim(Ca_sim);
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stim_delay(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_release, -pre_stim(2)*p_release];
        Fused(i) = pre_stim(2)*p_release;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_docking_transient,k_undocking_transient,syt3,syt7,Ca_sim/Ca_spike), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end

    state = [SS; state];
    ts = [delta_t; ts];
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


function dydt = dSS(~,state,k_docking,k_undocking,k_docking_transient,k_undocking_transient,syt3,syt7,Ca_res)

    
    dydt(1) = -state(1)*(k_docking + k_docking_transient*Ca_res) + state(2)*(k_undocking + k_undocking_transient*Ca_res/5);
    dydt(2) = state(1)*(k_docking + k_docking_transient*Ca_res) - state(2)*(k_undocking + k_undocking_transient*Ca_res/5);
    dydt = transpose(dydt);

end


function dydt = dState(t,state,k_docking,k_undocking,k_docking_transient,k_undocking_transient,syt3,syt7,Ca_sim_norm)
    
    Ca = Ca_sim_norm(round(t)+1);

    dydt(1) = -state(1)*(k_docking + k_docking_transient*Ca) + state(2)*(k_undocking + k_undocking_transient*Ca/5);
    dydt(2) = state(1)*(k_docking + k_docking_transient*Ca) - state(2)*(k_undocking + k_undocking_transient*Ca/5);
    dydt = transpose(dydt);

end