Ca_rest = 5e-8;
Ca_spike = 1.5e-5;
Ca_residual = 250e-9;
T_Ca_decay = 40; %ms
 
delta_t = 1e-2;
%stimulus_times = linspace(0,500*10,10);
%stimulus_times = [0 50 100 150 200 250 300 350 400 450];
%stimulus_times = [0 20 40 60 80 100 120 140 160 180];
stimulus_times = [0 10];
%stimulus_times = [0];
max_time = stimulus_times(end)+300;

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

ts = linspace(0,max_time,max_time/delta_t + 1);

k_on_1 = 1.4e5; %M^-1ms^-1 Kobbersmed
k_on_3 = 3e5; %Hui
%k_on_3 = 0; %syt3 KO
k_on_7 = 7.333e3; %Knight
%k_on_7 = 4.1e4;
%k_on_7 = 0; %syt7 KO
k_off_1 = 4; %Kobbersmed
k_off_3 = 1.5; %kobbersmed/Sugita
k_off_7 = 1.1e-2;
%k_off_7 = 6.15e-2; %Kobbersmed 
L_plus = 3.5e-7; %ms^-1
O = 27.978;
T = 1; %no effect of syt3 on release
%S = 1311; %est from Arrhenius
%S = 510.26; %best fit from Kobbersmed
S = 1; %no effect of syt7 on release
b1 = .5;
b3 = .5;
b7 = .5;
k_refill = 0.001;
%CDR = 30;
CDR = 1;
k_unprime = .2;
%K_A_prime = 1.7e-6; %Knight Ca1/2 syt7
K_A_prime = 55e-9; %Kobbersmed fits
k_prime = 0.002;
K_A_unprime = 7e-6; %Hui Ca1/2 syt3

syts = 3; %choose which syts you want to use, valid inputs [3, 7, 37]
new_params = 1; %set to one when testing new parameters to enable calculation of steady state

if new_params == 1
    
    filename = 'Three_sensor.xlsx'; %file containing rate and calcium dependence matricies
    k = [k_on_1; k_on_3; k_on_7; k_off_1; k_off_3; k_off_7; L_plus; O; T; S; b1; b3; b7; k_refill; CDR]; 
    writematrix(k,filename,'Sheet',1,'Range','B68'); %pass parameter values for calculation of the rate matrix
    
%     writematrix(k_refill,filename,'Range','B81');
%     writematrix(L_plus,filename,'Range','B74');
    
    Ca_ind = readmatrix(filename,'Sheet',3,'Range','B2:I9','UseExcel',1); %useexcel must be true for calculation to be made within the spreadsheet
    Ca_dep = readmatrix(filename,'Sheet',3,'Range','L2:S9','UseExcel',1);

    state_0 = zeros(8,1);
    state_0(2) = 1; %start with all sites in empty state

    %create rate matrix for steady state determination
    rate_matrix = Ca_rest*Ca_dep + Ca_ind;
    
    r = 1-Ca_rest^2/(Ca_rest^2 + K_A_prime^2);
    rate_matrix(2,3) = r*k_unprime;
    u = 1-Ca_rest^2/(Ca_rest^2 + K_A_unprime^2);
    rate_matrix(3,2) = u*k_prime;
    
    rate_matrix = rate_matrix - diag(sum(rate_matrix));

    %solution to ODEs for steady state determination, this should always be run
    %first when testing a new set of parameters
    [t0,state] = ode23s(@(t,state) SSvectorized(t,state,rate_matrix), ts, state_0);
    SS = state(end,:);
    SS(1) = 0;
    save('SS.mat','SS');
    
    %writematrix(L_plus,filename,'Sheet',1,'Range','B74');
    
    Ca_ind = readmatrix(filename,'Sheet',3,'Range','B2:I9','UseExcel',1); %useexcel must be true for calculation to be made within the spreadsheet
    Ca_dep = readmatrix(filename,'Sheet',3,'Range','L2:S9','UseExcel',1);
    
    save('SS.mat', 'SS');
    save('Ca_ind.mat', 'Ca_ind');
    save('Ca_dep.mat', 'Ca_dep');
    
else
    
    SS = matfile('SS.mat').SS;
    Ca_ind = matfile('Ca_ind.mat').Ca_ind;
    Ca_dep = matfile('Ca_dep.mat').Ca_dep;

end

Ca_sim = zeros(1,length(ts));
Ca_sim = Ca_sim + Ca_rest;
    
for t = 1:length(stimulus_times) %simulate calcium influx
   
    spike_start_index = round(stimulus_times(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one
    spike_peak_index = round((stimulus_times(t)+mu)/delta_t) + 1;
    
    Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
    
    Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

end

%solution to ODEs for stimuli

syt_selections = [3 7 37];
        
if sum(syt_selections == syts) ~= 1
    select = menu('Which syts would you like to use?',3,7,37);
    syts = syt_selections(select);
end

switch syts
    
    case 37
        %syt3 and 7
        [t,sol] = ode15s(@(t,sol) vectorized37(t,sol,Ca_sim,Ca_dep,Ca_ind,k_unprime,K_A_prime,k_prime,K_A_unprime,delta_t), ts, SS);
    
    case 3
        %syt3
        [t,sol] = ode15s(@(t,sol) vectorized3(t,sol,Ca_sim,Ca_dep,Ca_ind,k_unprime,k_prime,K_A_unprime,delta_t), ts, SS);
    
    case 7
        %syt7
        [t,sol] = ode15s(@(t,sol) vectorized7(t,sol,Ca_sim,Ca_dep,Ca_ind,k_unprime,K_A_prime,k_prime,delta_t), ts, SS);
    
    otherwise
        disp('ERROR')
        
end

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

function dydt = vectorized37(t,state,Ca_sim,Ca_dep,Ca_ind,k_unprime,K_A_prime,k_prime,K_A_unprime,delta_t)

    Ca = Ca_sim(round(t/delta_t)+1); 

    rate_matrix = Ca*Ca_dep + Ca_ind;
    rate_matrix(2,3) = (1-Ca^2/(Ca^2 + K_A_prime^2))*k_unprime; %inhibition of unpriming
    rate_matrix(3,2) = (1+Ca^2/(Ca^2 + K_A_unprime^2))*k_prime; %acceleration of priming
    rate_matrix = rate_matrix - diag(sum(rate_matrix));
    
    dydt = rate_matrix*state; %runs slower if rate_matrix is set to sparse(rate_matrix)
    dydt(2) = dydt(2) + dydt(1);
end

function dydt = vectorized3(t,state,Ca_sim,Ca_dep,Ca_ind,k_unprime,k_prime,K_A_unprime,delta_t)

    Ca = Ca_sim(round(t/delta_t)+1); 

    rate_matrix = Ca*Ca_dep + Ca_ind;
    rate_matrix(2,3) = k_unprime; %no syt7 effect
    rate_matrix(3,2) = (1+Ca^2/(Ca^2 + K_A_unprime^2))*k_prime;
    rate_matrix = rate_matrix - diag(sum(rate_matrix));
    
    dydt = rate_matrix*state; %runs slower if rate_matrix is set to sparse(rate_matrix)
    dydt(2) = dydt(2) + dydt(1);
end

function dydt = vectorized7(t,state,Ca_sim,Ca_dep,Ca_ind,k_unprime,K_A_prime,k_prime,delta_t)

    Ca = Ca_sim(round(t/delta_t)+1); 

    rate_matrix = Ca*Ca_dep + Ca_ind;
    rate_matrix(2,3) = (1-Ca^2/(Ca^2 + K_A_prime^2))*k_unprime; 
    rate_matrix(3,2) = k_prime; %no syt3 effect
    rate_matrix = rate_matrix - diag(sum(rate_matrix));
    
    dydt = rate_matrix*state; %runs slower if rate_matrix is set to sparse(rate_matrix)
    dydt(2) = dydt(2) + dydt(1);
end