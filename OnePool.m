global p_release
global k_docking
global k_undocking
global delta_t
global SS

p_release = .8;
delta_t = 0.01;
k_docking = .0012;
k_undocking = .0004;

state_0 = [1,0]; %empty,filled
t_SS = 10000;
[~,state] = ode45(@(t,state) dState(t,state,k_docking,k_undocking), [0 t_SS], state_0);

SS = state(end,:);

%CF data
train_recovery_data_WT = [0.214594108 0.384408264 0.638176057 0.817518904 0.887058714 0.872159256 0.964188913 0.963314927 0.982736016 0.981537505]; %set 8/23/20
train_recovery_data_SEMs_WT = [0.035160387 0.040212604 0.043887001 0.027611691 0.025484558 0.05457874 0.03331073 0.024722031 0.011843753 0.011107924]; %set 8/23/20

% train_recovery_data_WT_Ca_block = [.133 .301 .535 .706 .853 .918 .950 1.03 1.01 1.03];
% train_recovery_data_SEMs_WT_Ca_block = [0.011220544	0.02984717 0.052266679 0.048872772 0.029092336 0.026307037 0.021643624 0.012835957 0.009772508 0.027547534];

% train_recovery_data_syt3KO_Ca_block = [0.095687504 0.1867 0.483622448 0.832322146 1.02638655 1.15260018	1.049615847	0.890841363	1.010482505	1.002926784];
% train_recovery_data_SEMs_syt3KO_Ca_block = [0.020921145	0.025491913	0.053612745	0.036642093	0.023708192	0.058597074	0.020695876	0.066537454	0.012680184	0.017007916];
%
train_recovery_data_syt3KO = [0.115660845 0.165504058 0.336442556 0.610561214 0.750744064 0.896371099 0.952219223 1.043590101 1.040256315 1.009667422];
train_recovery_data_SEMs_syt3KO = [0.020201196 0.032159398 0.042618698 0.050037539 0.031454932 0.021833684 0.010196577 0.023865478 0.020565531 0.015465173];

% train_recovery_data = train_recovery_data_syt3KO_Ca_block;
% train_recovery_data_SEMs = train_recovery_data_SEMs_syt3KO_Ca_block;

train_recovery_data = train_recovery_data_WT;
train_recovery_data_SEMs = train_recovery_data_SEMs_WT;

% train_recovery_data = train_recovery_data_WT_Ca_block;
% train_recovery_data_SEMs = train_recovery_data_SEMs_WT_Ca_block;

% train_recovery_data = train_recovery_data_syt3KO;
% train_recovery_data_SEMs = train_recovery_data_SEMs_syt3KO;


% trainSS_data_WT = [0.8162 0.4297 0.2587 0.00896]; %set 8/23/20
% trainSS_data_SEMs_WT = [0.0048 0.0037 0.0077 0.0027]; %set 8/23/20

% trainSS_data_WT_Ca_block = [.722 .241 .104 .0138]; 

trainSS_data_syt3KO = [0.7419  0.2852 0.1350 0.0676]; %set 8/23/20
trainSS_data_SEMs_syt3KO = [0.0021 0.0027 0.0022 0.0010]; %set 8/23/20

% trainSS_data_syt3KO_Ca_block = [0.761628403 0.276455062 0.127877809 0.016273454];
% 
% trainSS_data = trainSS_data_syt3KO_Ca_block;

% trainSS_data = trainSS_data_WT;
% trainSS_data_SEMs = trainSS_data_SEMs_WT;

% trainSS_data = trainSS_data_WT_Ca_block;

trainSS_data = trainSS_data_syt3KO;
trainSS_data_SEMs = trainSS_data_SEMs_syt3KO;


% 
% PPR_data_WT = [0.22 0.32 0.39 0.42 0.50 0.58 0.70 0.75 0.80 0.84 0.87 0.94 0.97 0.98 1.00 1.01]; %set 8/23/20
% PPR_data_SEMs_WT = [0.02 0.03 0.03 0.03 0.02 0.02 0.03 0.03 0.02 0.02 0.02 0.01 0.01 0.02 0.01 0.01]; %set 8/23/20

% PPR_data_WT_Ca_block = [.103 .194 .272 .307 .346 .365 .417 .477 .503 .620 .615 .850 .936 .965 .980 1.04]; 
% PPR_data_SEMs_WT_Ca_block = [0.023707394 0.031631686 0.036341527 0.024351258 0.020447046 0.022314029 0.02242531	0.024670909	0.024094207	0.020289501	0.019038797	0.012231847	0.012046944	0.006259344	0.010283128	0.00864487];

PPR_data_syt3KO = [0.15 0.21 0.27 0.30 0.34 0.40 0.55 0.68 0.74 0.79 0.83 0.90 0.96 0.98 0.99 1.03]; %set 8/23/20
PPR_data_SEMs_syt3KO = [0.02 0.02 0.03 0.04 0.04 0.05 0.05 0.02 0.03 0.03 0.03 0.02 0.01 0.01 0.01 0.01]; %set 8/23/20

% PPR_data_syt3KO_Ca_block = [0.151147669	0.184327088	0.238073419	0.278179684	0.349569587	0.371625114	0.4931594 0.601356914 0.617759941 0.676810496 0.743110597 0.876359759 0.944711679 0.973497529 0.976220011 1.017493162];
% PPR_data_SEMs_syt3KO_Ca_block = [0.026678818 0.023828808 0.033945846 0.031585123 0.030657408 0.032225669 0.022723458 0.051859004 0.053610984 0.051481764 0.040548383 0.028425046 0.021286098 0.018492625 0.01164004	0.010175575];
% 
% PPR_data = PPR_data_syt3KO_Ca_block;
% PPR_data_SEMs = PPR_data_SEMs_syt3KO_Ca_block;

% PPR_data = PPR_data_WT;
% PPR_data_SEMs = PPR_data_SEMs_WT;

% PPR_data = PPR_data_WT_Ca_block;
% PPR_data_SEMs = PPR_data_SEMs_WT_Ca_block;

PPR_data = PPR_data_syt3KO;
PPR_data_SEMs = PPR_data_SEMs_syt3KO;


%begin calls
%PPR sim
ISI_PPR = [10 20 35 50 75 100 200 350 500 750 1000 2000 3500 5000 7500 10000];
PPR = zeros(1,length(ISI_PPR));

for j = 1:length(ISI_PPR)
    
    stimulus_times = [0 ISI_PPR(j)];

    max_time = stimulus_times(end) + stimulus_times(2)*3;

    [~, ~, Fused] = stim_sim(stimulus_times, max_time);
    Fused = Fused/Fused(1);
    PPR(j) = Fused(2);

end

%train recovery sim
stimulus_times = linspace(0,20*19,20); %20 stims @50hz

max_time = stimulus_times(end) + stimulus_times(2)*10;

[ts, state, Fused] = stim_sim(stimulus_times, max_time);

post_train = state(ts==stimulus_times(end),:);

ISI_train_recovery = [50 100 200 350 500 750 1000 2000 5000 10000];
train_recovery = zeros(1,length(ISI_train_recovery));

for j = 1:length(ISI_train_recovery)

    [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking), [0 ISI_train_recovery(j)], post_train);
    train_recovery(j) = out(end,2)*p_release/Fused(1);

end

%SS/freq sim

freq = [1 10 20 50];
trainSS = zeros(1,length(freq));

for j = 1:length(freq)
    
    stimulus_times = linspace(0, 1000/freq(j)*19, 20);

    max_time = stimulus_times(end) + stimulus_times(2)*10;

    [~, ~, Fused] = stim_sim(stimulus_times, max_time);

    Fused = Fused/Fused(1);
    trainSS(j) = mean(Fused(end-5:end)); %take last 6 values to get avg SS

end

%plots

subplot(3,1,1)
semilogx(ISI_PPR,PPR,'-b.')
title('PPR simulation and data')
xlabel('Time(ms)')
ylabel('PPR')
hold on
errorbar(ISI_PPR,PPR_data,PPR_data_SEMs,'-r')
legend('Simulation','Data','Location','southeast')
hold off

subplot(3,1,2)
semilogx(ISI_train_recovery,train_recovery,'-b.')
title('Train recovery simulation and data')
xlabel('Time(ms)')
ylabel('Recovery')
hold on
errorbar(ISI_train_recovery,train_recovery_data,train_recovery_data_SEMs,'-r')
legend('Simulation','Data','Location','southeast')
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

%begin function definitions
function [ts, state, Fused] = stim_sim(stimulus_times, max_time)
    
    global p_release
    global k_docking
    global k_undocking
    global delta_t
    global SS
    
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stim_delay(end)];
    ts = 0;
    Fused = zeros(length(stimulus_times),1);
    
    state = SS;
    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        released = pre_stim(2)*p_release;
        post_stim = pre_stim + [released, -released];
        Fused(i) = released;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end
    
    state = [SS; state];
    ts = [delta_t; ts];
    
end


function dydt = dState(~,state,k_docking,k_undocking)

    dydt(1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2) = state(1)*k_docking - state(2)*k_undocking;
    dydt = transpose(dydt);

end