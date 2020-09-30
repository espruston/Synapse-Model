%stimulus_times = [0 10];
%stimulus_times = linspace(0,20*19,20); %20 stims @50hz

K_D_3 = 3e5; %M^-1ms^-1 Hui
K_D_7 = 7.333e3; %M^-1ms^-1 Brandt/Knight 

% hold values, PPR R^2 = .9899
% k_docking = .0004;
% k_undocking = .0001;
% k_maturation = .015;
% k_dematuration = .01;

k_docking = .0004;
k_undocking = .0001;
k_maturation = .015;
k_dematuration = .01;

CDR = 0;
Facil = 0;

C_3 = 0;
C_7 = 0;

p_immature = .1;
p_mature = 1;

Ca_rest = 5e-8; %M
Ca_spike = 2e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms

delta_t = 1e-2; %ms

t_SS = 10000; %ms
ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [0; 1; 0]; %start all vesicles in immature docked state 

[t0,SS] = ode15s(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration,Ca_rest), [0 t_SS], state_0);

freq = [1 10 20 50];
train_SS = zeros(1,length(freq));

for j = 1:length(freq)
    state = SS(end,:);
    stimulus_times = linspace(0, 1000/freq(j)*19, 20);
    
    max_time = stimulus_times(end) + stimulus_times(2)*10;
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stim_delay(end)];

    ts = 0;
    Fused = zeros(length(stimulus_times),1);

    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_immature + pre_stim(3)*p_mature -pre_stim(2)*p_immature -pre_stim(3)*p_mature];
        Fused(i) = pre_stim(2)*p_immature + pre_stim(3)*p_mature;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,Ca_rest), [0 stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1); t+ts(end)];
    end

    Fused = Fused/Fused(1);
    train_SS(j) = mean(Fused(end-5:end));

end

WT_Ca_block = [.722 .241 .104 .0138]; 

R = corrcoef(train_SS, WT_Ca_block);
disp(R.^2)

plot(freq, train_SS, 'bx')
title('Steady state of train simulation and data')
xlabel('Time (ms)')
ylabel('Train steady state')
hold on
plot(freq, WT_Ca_block, 'go')
legend('Simulation','Data')
hold off

function dydt = dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration, Ca_rest)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration, Ca_rest)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
    
end