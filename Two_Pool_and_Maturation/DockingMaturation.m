%stimulus_times = [0 10];
stimulus_times = linspace(0,1000*19,20); %20 stims
stimulus_times = linspace(0,50*19,20); %20 stims

K_D_3 = 3e5; %M^-1ms^-1 Hui
K_D_7 = 7.333e3; %M^-1ms^-1 Brandt/Knight 

% hold values, PPR R^2 = .9724
% k_docking = .03;
% k_undocking = .0001;
% k_maturation = .000965;
% k_dematuration = .0001;

k_docking = .05;
k_undocking = 0;
k_maturation = .002;
k_dematuration = 0;

CDR = 0;
Facil = 0;

C_3 = 0;
C_7 = 0;

p_immature = .5;
p_mature = 1;

Ca_rest = 5e-8; %M
Ca_spike = 2e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms

delta_t = 1e-2; %ms

t_SS = 10000000; %ms
ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

state_0 = [0; 1; 0]; %start all vesicles in immature docked state 

[t0,SS] = ode15s(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration), [0 t_SS], state_0);

state = SS(end,:);

max_time = stimulus_times(end) + stimulus_times(2)*10;
stim_delay = diff(stimulus_times);
stim_delay = [stim_delay max_time-stim_delay(end)];

ts = 0;
Fused = zeros(length(stimulus_times),1);

for i = 1:length(stim_delay)
    
    pre_stim = state(end,:);
    post_stim = pre_stim + [pre_stim(2)*p_immature + pre_stim(3)*p_mature -pre_stim(2)*p_immature -pre_stim(3)*p_mature];
    Fused(i) = pre_stim(2)*p_immature + pre_stim(3)*p_mature;
    [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration), [0 stim_delay(i)], post_stim);
    
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

subplot(2,2,1)
plot(Fused,'b.','MarkerSize', 20)
xlabel('Pulse #')
ylabel('Peak EPSC')
set(gca,'xlim',[1 length(Fused)])
set(gca,'ylim',[0 1])
hold on

hz_1 = [1 0.904585375 0.906090801 0.885396962 0.888286409 0.882153581 0.868080327 0.870438248 0.858248517 0.854839359 0.840352061 0.845592803 0.847296096 0.844539428 0.829714224 0.830002447 0.812125568 0.815760427 0.81009733 0.799518304];
hz_10 = [1 0.616949107 0.628326956 0.604464053 0.565850677 0.537103779 0.525309915 0.503692205 0.474574723 0.470769015 0.461639789 0.45151103 0.450073867 0.451884135 0.44312584 0.429802182 0.433756546 0.431020365 0.423930097 0.416302277];
hz_20 = [1 0.489671267 0.444366914 0.457655745 0.430402833 0.416929064 0.38834963 0.351531591 0.348392421 0.329990273 0.31629381 0.305756525 0.3022066 0.288966797 0.283932027 0.266221279 0.276873203 0.265228816 0.262831634 0.258863279];
hz_50 = [1 0.344586239 0.28753639 0.241585061 0.22098407 0.212199666 0.198432318 0.177849609 0.158522769 0.144448699 0.131731778 0.119584277 0.111943939 0.106093112 0.098933385 0.095593041 0.08975884 0.087444448 0.082325078 0.083824248];
 
    
plot(hz_1)
plot(hz_10)
plot(hz_20)
plot(hz_50)

hold off

subplot(2,2,2)
plot(ts, state(:,1),'-b')
xlabel('Time (ms)')
ylabel('Empty sites')
set(gca,'xlim',[ts(1) max_time])
set(gca,'ylim',[0 1])

subplot(2,2,3)
plot(ts, state(:,2),'-b')
xlabel('Time (ms)')
ylabel('Immature Vesicles')
set(gca,'xlim',[ts(1) max_time])
set(gca,'ylim',[0 1])

subplot(2,2,4)
plot(ts, state(:,3),'-b')
xlabel('Time (ms)')
ylabel('Mature Vesicles')
set(gca,'xlim',[ts(1) max_time])
set(gca,'ylim',[0 1])

function dydt = dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end