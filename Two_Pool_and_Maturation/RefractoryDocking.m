%stimulus_times = [0 10];
stimulus_times = linspace(0,20*19,20); %20 stims @50hz

K_D_3 = 3e5; %M^-1ms^-1 Hui
K_D_7 = 7.333e3; %M^-1ms^-1 Brandt/Knight 

% hold values, PPR R^2 = .9724
% k_docking = .03;
% k_undocking = .0001;
% k_maturation = .000965;
% k_dematuration = .0001;

k_refractory = 0;
k_docking = .03;
k_undocking = .001;
k_maturation = 0;
k_dematuration = 0;

CDR = 0;
Facil = 0;

C_3 = 0;
C_7 = 0;

p_immature = .9;
p_mature = .9;

Ca_rest = 5e-8; %M
Ca_spike = 2e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms

delta_t = 1e-2; %ms

t_SS = 10000; %ms
ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [0; 0; 1; 0]; %start all vesicles in immature docked state 

[t0,SS] = ode15s(@(t,state) dSS(t,state,k_refractory,k_docking,k_undocking,k_maturation,k_dematuration,Ca_rest), [0 t_SS], state_0);

state = SS(end,:);

max_time = stimulus_times(end) + stimulus_times(2)*10;
stim_delay = diff(stimulus_times);
stim_delay = [stim_delay max_time-stim_delay(end)];

ts = 0;
Fused = zeros(length(stimulus_times),1);

for i = 1:length(stim_delay)
    
    pre_stim = state(end,:);
    post_stim = pre_stim + [pre_stim(3)*p_immature + pre_stim(4)*p_mature, 0, -pre_stim(3)*p_immature, -pre_stim(4)*p_mature];
    Fused(i) = pre_stim(3)*p_immature + pre_stim(4)*p_mature;
    [t,out] = ode45(@(t,state) dState(t,state,k_refractory,k_docking,k_undocking,k_maturation,k_dematuration,Ca_rest), [0 stim_delay(i)], post_stim);
    
    state = [state(1:end-1,:); out];
    
    ts = [ts(1:end-1); t+ts(end)];
end

Fused = Fused/Fused(1);


y = exp(-1*ts/.1)-exp(-1*ts/2);
y = -y/min(y);
EPSC_sim = zeros(length(ts),1);

for i = 1:length(stimulus_times)
    start_index = find(ts==stimulus_times(i));
    EPSC_sim(start_index:end) = EPSC_sim(start_index:end) + y(1:end-start_index+1)*Fused(i);
end

ts = [-10; -delta_t; ts];
state = [SS(end,:); SS(end,:); state];
%figure

subplot(3,2,1)
plot(ts, [0; 0; EPSC_sim])
xlabel('Time (ms)')
ylabel('EPSC (norm)')
set(gca,'xlim',[ts(1) max_time])

subplot(3,2,2)
plot(ts, state(:,1))
xlabel('Time (ms)')
ylabel('Refractory sites')
set(gca,'xlim',[ts(1) max_time])

subplot(3,2,3)
plot(ts, state(:,2))
xlabel('Time (ms)')
ylabel('Empty available sites')
set(gca,'xlim',[ts(1) max_time])

subplot(3,2,4)
plot(ts, state(:,3))
xlabel('Time (ms)')
ylabel('Docked sites')
set(gca,'xlim',[ts(1) max_time])

subplot(3,2,5)
plot(ts, state(:,4))
xlabel('Time (ms)')
ylabel('Mature sites')
set(gca,'xlim',[ts(1) max_time])

function dydt = dSS(t,state,k_refractory,k_docking,k_undocking,k_maturation,k_dematuration, Ca_rest)
    
    dydt(1,1) = -state(1)*k_refractory;
    dydt(2,1) = state(1)*k_refractory - state(2)*k_docking + state(3)*k_undocking;
    dydt(3,1) = state(2)*k_docking - state(3)*k_undocking - state(3)*k_maturation + state(4)*k_dematuration;
    dydt(4,1) = state(3)*k_maturation - state(4)*k_dematuration;

end

function dydt = dState(t,state,k_refractory,k_docking,k_undocking,k_maturation,k_dematuration, Ca_rest)
    
    dydt(1,1) = -state(1)*k_refractory;
    dydt(2,1) = state(1)*k_refractory - state(2)*k_docking + state(3)*k_undocking;
    dydt(3,1) = state(2)*k_docking - state(3)*k_undocking - state(3)*k_maturation + state(4)*k_dematuration;
    dydt(4,1) = state(3)*k_maturation - state(4)*k_dematuration;

end