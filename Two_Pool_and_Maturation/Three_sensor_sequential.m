Ca_rest = 5e-8; %M, Skyler poster
Ca_spike = 5e-5; %M
Ca_residual = 250e-9; %M
T_Ca_decay = 40; %ms
 
delta_t = 1e-2;
%stimulus_times = linspace(0,500*10,10);
%stimulus_times = [0 50 100 150 200 250 300 350 400 450];
%stimulus_times = [0 20 40 60 80 100 120 140 160 180];
%stimulus_times = [10 20 30 40 50 60 70 80 90 100];
stimulus_times = [0 10];
%stimulus_times = [0];
max_time = stimulus_times(end)+500;

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

ts = linspace(0,max_time,max_time/delta_t + 1);

k_on_1 = 1.4e5; %M^-1ms^-1 Ca binding rate, syt1, Kobbersmed
k_on_3 = 3e5; %M^-1ms^-1 membrane binding rate, syt3, Hui
%k_on_3 = 0; %turn off syt3
k_on_7 = 7.333e3; %M^-1ms^-1 membrane binding rate, syt7, Knight
%k_on_7 = 0; %tun off syt7

k_off_1 = 4; %ms^-1 Ca unbinding rate, Kobbersmed
% k_off_3 = 1.28e-2; %ms^-1 membrane unbinding rate, Hui
% k_off_7 = 1.1e-2; %ms^-1 membrane unbinding rate, Knight\
b = .5;

L_plus = 3.5e-7; %ms^-1 Kobbersmed
%writematrix(L_plus,filename,'Sheet',1,'Range','B74');

O = 27.978; %fusion rate constant for syt1, Kobbersmed

k_basal = 2e-3; %ms^-1
%k_prime = k_on_3;
%k_prime = 5e3;
k_prime = 0;
k_unprime = 5e-2; %ms^-1
k_fill = 8e3;
%k_unfill_inhib = k_on_7;
k_unfill_inhib = 0;
k_unfill = .12;


vars = [L_plus; O; k_on_1; k_off_1; b; k_basal; k_prime; k_unprime; k_fill; k_unfill_inhib; k_unfill];


Ca_sim = zeros(1,length(ts));
Ca_sim = Ca_sim + Ca_rest;
    
for t = 1:length(stimulus_times) %simulate calcium influx
   
    spike_start_index = round(stimulus_times(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one
    spike_peak_index = round((stimulus_times(t)+mu)/delta_t) + 1;
    
    Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
    
    Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

end

filename = 'Sequential.xlsx'; %file containing rate and calcium dependence matricies
writematrix(vars,filename,'Sheet',1,'Range','B13'); %pass parameter values for calculation of the rate matrix

Ca_ind = readmatrix(filename,'Sheet',1,'Range','B2:J10','UseExcel',1); %useexcel must be true for calculation to be made within the spreadsheet
Ca_dep = readmatrix(filename,'Sheet',1,'Range','M2:U10','UseExcel',1);

state_0 = zeros(9,1);
state_0(2) = 1; %start with all vesicles in reserve

%create rate matrix for steady state determination
rate_matrix = Ca_rest*Ca_dep + Ca_ind;
rate_matrix = rate_matrix - diag(sum(rate_matrix));

%solution to ODEs for steady state determination, this should always be run
%first when testing a new set of parameters
[t0,state] = ode23s(@(t,state) SSvectorized(t,state,rate_matrix), ts, state_0);

SS = state(end,:);
SS(1) = 0;

[t,sol] = ode15s(@(t,sol) vectorized(t,sol,Ca_sim,Ca_dep,Ca_ind,delta_t), ts, SS);

Fused = sol(:,1);
y = exp(-1*ts/.1)-exp(-1*ts/2);
dFused = [0; diff(Fused)];
EPSC_sim = zeros(1,length(dFused));

for i = 1:length(dFused)
EPSC_sim(i:end) = EPSC_sim(i:end) + y(1:end-i+1)*dFused(i);
end
EPSC_sim = -1*EPSC_sim/min(EPSC_sim(1:(stimulus_times(1)+2)/delta_t)); %normalize to the first EPSC under the assumption that the peak occurs within 2ms of the stimulus

plot([-10 ts], [0 EPSC_sim])
xlim([-10 max_time])


function dydt = SSvectorized(t,state,rate_matrix)
    dydt = rate_matrix*state;
    dydt(2) = dydt(2) + dydt(1);
end

function dydt = vectorized(t,state,Ca_sim,Ca_dep,Ca_ind,delta_t)

    Ca = Ca_sim(round(t/delta_t)+1); 

    rate_matrix = Ca*Ca_dep + Ca_ind;
    rate_matrix = rate_matrix - diag(sum(rate_matrix));
    
    dydt = rate_matrix*state; %runs slower if rate_matrix is set to sparse(rate_matrix)
    dydt(2) = dydt(2) + dydt(1);
end