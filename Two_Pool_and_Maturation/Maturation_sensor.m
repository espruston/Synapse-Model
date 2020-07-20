%maturation_v2 assumes a maturation model where the maturation can be
%accelerated by calcium from a pulse

k_refill = 0; %vesicles ms-1

k_mature_basal = 2; % basal rate of vesicle maturation, ms-1
k_mature = 3e5; % syt3, M-1ms-1, Hui
k_unmature = .378; %ms-1, Hui

k_prime = 7.33e3; % syt7, M-1ms-1, Knight
k_unprime = .011; %ms-1

k_on = 1.4e5; %ca binding rate, syt1, M-1ms-1
k_off = 4; %ca unbinding rate, syt1, ms-1 

b = 0.5; %cooperativity of Ca binding syt1

f = 28;

Ca_rest = 5e-8;
Ca_spike = 2.5e-5;
Ca_residual = 250e-9;
T_Ca_decay = 40; %ms
 
delta_t = 1e-2;
max_time = 500;
stimulus_times = [0 20 40 60 80 100 120 140 160 180];
%stimulus_times = [0];

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

ts = linspace(0,max_time,max_time/delta_t + 1);
I = eye(21);

new_params = 1; %set to one when testing new parameters to enable calculation of steady state

if new_params == 1
    %values for simulation of steady state
    k_fuse_basal = 0;
    k_refill = 0;
    M_plus = 0;
    P_plus = 0;

    filename = 'rate_matrix_Maturation_sensor.xlsx'; %file containing rate and calcium dependence matricies
    k = [k_refill;k_mature_basal;k_mature;k_unmature;k_prime;k_unprime;k_on;k_off;b;f;M_plus;P_plus]; 
    writematrix(k,filename,'Sheet',2,'Range','B48'); %pass parameter values for calculation of the rate matrix

    Ca_independent_matrix = readmatrix(filename,'Sheet',2,'Range','B2:V22','UseExcel',1); %useexcel must be true for calculation to be made within the spreadsheet
    Ca_dependence_matrix = readmatrix(filename,'Sheet',2,'Range','B25:V45','UseExcel',1);
    
%     %these matricies are sparse due to relatively low connectivity between states
%     Ca_independent_matrix = triu(M); %upper triangular by design
%     Ca_dependence_matrix = tril(M); %should be lower triangular by design

    state_0 = zeros(21,1);
    state_0(2) = 1; %start all vesicles in immature state

    %create rate matrix for steady state determination
    rate_matrix = Ca_rest*Ca_dependence_matrix + Ca_independent_matrix;
    rate_matrix = rate_matrix - sum(rate_matrix).*I;

    %solution to ODEs for steady state determination, this should always be run
    %first when testing a new set of parameters
    [t0,state] = ode45(@(t,state) SSvectorized(t,state,rate_matrix), ts, state_0);

    %define each 'type' of vesicle for SS determination, summation occurs along axis 2
    Immature = state(:,2);
    Mature = sum(state(:,[3 7:11]),2); %mature vesicle count should include those which are Ca bound
    Primed = sum(state(:,[4 12:16]),2);
    Superprimed = sum(state(:,[6 17:21]),2); %state(:,5) contains reserve
    Fused = sum(state(:,1),2);

    SS = zeros(21,1);

    SS([1 2 3 4 6]) = [Fused(end) Immature(end) Mature(end) Primed(end) Superprimed(end)];

    save('SS.mat','SS');
    %turn on fusion, simulate calcium influx
    %values for simulation of stims
    k_fuse_basal = 3.5e-4;
    %k_refill = 0;
    M_plus = 3.5e-4;
    P_plus = 3.5e-4;
    
    %update rate matrix for non SS solution
    %writematrix(k_refill,filename,'Sheet',2,'Range','B48');
    k = [M_plus;P_plus]; 
    writematrix(k,filename,'Sheet',2,'Range','B58'); %pass parameter values for calculation of the rate matrix

    Ca_independent_matrix = readmatrix(filename,'Sheet',2,'Range','B2:V22','UseExcel',1); %useexcel must be true for calculation to be made within the spreadsheet
    Ca_dependence_matrix = readmatrix(filename,'Sheet',2,'Range','B25:V45','UseExcel',1);
    
    save('SS.mat', 'SS');
    save('Ca_independent.mat', 'Ca_independent_matrix');
    save('Ca_dependence.mat', 'Ca_dependence_matrix');
    
else
    
    SS = matfile('SS.mat').SS;
    Ca_independent_matrix = matfile('Ca_independent.mat').Ca_independent_matrix;
    Ca_dependence_matrix = matfile('Ca_dependence.mat').Ca_dependence_matrix;

end

Ca_sim = zeros(1,length(ts));
Ca_sim = Ca_sim + Ca_rest;
    
for t = 1:length(stimulus_times) %simulate calcium influx
   
    spike_start_index = stimulus_times(t)/delta_t + 1; %if 1st stim is at t=0 index should be one
    spike_peak_index = (stimulus_times(t)+mu)/delta_t + 1;
    
    Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one
    
    Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

end

%solution to ODEs for stimuli
[t,sol] = ode45(@(t,sol) vectorized(t,sol,Ca_sim,Ca_dependence_matrix,Ca_independent_matrix,I,delta_t), ts, SS);

%define each 'type' of vesicle for plotting, summation occurs along axis 2
Immature = sol(:,2);
Mature = sum(sol(:,[3 7:11]),2); %mature vesicle count should include those which are Ca bound
Primed = sum(sol(:,[4 12:16]),2);
Superprimed = sum(sol(:,[6 17:21]),2); %state(:,5) contains reserve
Fused = sum(sol(:,1),2);

%calculate simulated EPSC as a double exponential
y = exp(-1*ts/.1)-exp(-1*ts/2);
dFused = [0; diff(Fused)];
EPSC_sim = zeros(1,length(dFused));

for i = 1:length(dFused)
EPSC_sim(i:end) = EPSC_sim(i:end) + y(1:end-i+1)*dFused(i);
end
EPSC_sim = -1*EPSC_sim/min(EPSC_sim(1:(stimulus_times(1)+2)/delta_t)); %normalize to the first EPSC under the assumption that the peak occurs within 2ms of the stimulus

subplot(4,2,1);
plot([-10 ts], [0 EPSC_sim]);
title('Simulated EPSC')
xlabel('time (ms)')
ylabel('EPSC_{normalized}')
xlim([-10 max_time])


subplot(4,2,2);
plot([-10 ts], [0; dFused/delta_t]);
title('dFused')
xlabel('time (ms)')
ylabel('dFused')
xlim([-10 max_time])

subplot(4,2,3);
semilogy([-10 ts], [Ca_rest Ca_sim]);
title('Intracellular Ca signal')
xlabel('time (ms)')
ylabel('Ca signal (M)')
xlim([-10 max_time])
ylim([4e-8, max(Ca_sim)*1.5])


subplot(4,2,4);
plot([-10 ts], [0; Fused]);
title('Fused')
xlabel('time (ms)')
ylabel('Fused Vesicles')
xlim([-10 max_time])

subplot(4,2,5);
plot([-10 ts], [SS(2); Immature]);
title('Immature')
xlabel('time (ms)')
ylabel('Immature Vesicles')
xlim([-10 max_time])
%ylim([0, max(EPSC)*1.2])

subplot(4,2,6);
plot([-10 ts], [SS(3); Mature]);
title('Mature')
xlabel('time (ms)')
ylabel('Mature Vesicles')
xlim([-10 max_time])
%ylim([0, max(EPSC)*1.2])

subplot(4,2,7);
plot([-10 ts], [SS(4); Primed]);
title('Primed')
xlabel('time (ms)')
ylabel('Primed Vesicles')
xlim([-10 max_time])

subplot(4,2,8);
plot([-10 ts], [SS(6); Superprimed]);
title('Superprimed')
xlabel('time (ms)')
ylabel('Superprimed Vesicles')
xlim([-10 max_time])

function dydt = SSvectorized(t,state,rate_matrix)
    dydt = rate_matrix*state;
end

function dydt = vectorized(t,state,Ca_sim,Ca_dependence_matrix,Ca_independent_matrix,I,delta_t)

    Ca = Ca_sim(round(t)/delta_t+1); 

    rate_matrix = Ca*Ca_dependence_matrix + Ca_independent_matrix;
    rate_matrix = rate_matrix - sum(rate_matrix).*I;
    
    dydt = rate_matrix*state; %runs slower if rate_matrix is set to sparse(rate_matrix)
end



% function partials = maturation_v2(t,y)
% 
%     Immature = y(1);
%     Mature = y(2);
%     Prime = y(3);
%     Superprime = y(4);
%     
%     Mature_1 = y(5);
%     Mature_2 = y(6);
%     Mature_3 = y(7);
%     Mature_4 = y(8);
%     Mature_5 = y(9);
%     
%     Prime_1 = y(10);
%     Prime_2 = y(11);
%     Prime_3 = y(12);
%     Prime_4 = y(13);
%     Prime_5 = y(14);
%     
%     Superprime_1 = y(15);
%     Superprime_2 = y(16);
%     Superprime_3 = y(17);
%     Superprime_4 = y(18);
%     Superprime_5 = y(19);
%     
%     Fused = y(20);
%     Ca = y(21);
%   
%     
%     
%     dImmature_ = k_refill*(1-Immature) - (Ca*k_mature + Ca*k_prime)*Immature + k_unmature*Mature + k_unprime*Prime;
%     dMature = Ca*k_mature + k_unprime*Superprime + k_off_1*Mature_1 - (k_immature + Ca*k_prime + 5*Ca*k_on_1 + M_plus)*Mature;
%     dPrime = Ca*k_prime*Immature + k_unprime*Superprime + k_off_1*Prime_1 - (k_unprime + Ca*k_mature + 5*Ca*k_on_1 + P_plus)*Prime;
%     dSuperprime = Ca*k_prime*Mature + Ca*k_mature*Prime + k_off_1*Superprime_1 - (k_unprime + k_unmature + 5*Ca*k_on_1 + (M_plus + P_plus))*Superprime;
%     
%     dMature_1 = 5*Ca*k_on_1*Mature + 2*b*k_off_1*Mature_2 - (k_off_1 + 4*Ca*k_on_1 + M_plus*f)*Mature_1;
%     dMature_2 = 4*Ca*k_on_1*Mature_1 + 3*b^2*k_off_1*Mature_3 - (2*b*k_off_1 + 3*Ca*k_on_1 + M_plus*f^2)*Mature_2;
%     dMature_3 = 3*Ca*k_on_1*Mature_2 + 4*b^3*k_off_1*Mature_4 - (3*b^2*k_off_1 + 2*Ca*k_on_1 + M_plus*f^3)*Mature_3;
%     dMature_4 = 2*Ca*k_on_1*Mature_3 + 5*b^4*k_off_1*Mature_5 - (4*b^3*k_off_1 + Ca*k_on_1 + M_plus*f^4)*Mature_4;
%     dMature_5 = Ca*k_on_1*Mature_4 - (5*b^3*k_off_1 + M_plus*f^5)*Mature_5;
%     
%     dPrime_1 = 5*Ca*k_on_1*Prime + 2*b*k_off_1*Prime_2 - (k_off_1 + 4*Ca*k_on_1 + P_plus*f)*Prime_1;
%     dPrime_2 = 4*Ca*k_on_1*Prime_1 + 3*b^2*k_off_1*Prime_3 - (2*b*k_off_1 + 3*Ca*k_on_1 + P_plus*f^2)*Prime_2;
%     dPrime_3 = 3*Ca*k_on_1*Prime_2 + 4*b^3*k_off_1*Prime_4 - (3*b^2*k_off_1 + 2*Ca*k_on_1 + P_plus*f^3)*Prime_3;
%     dPrime_4 = 2*Ca*k_on_1*Prime_3 + 5*b^4*k_off_1*Prime_5 - (4*b^3*k_off_1 + Ca*k_on_1 + P_plus*f^4)*Prime_4;
%     dPrime_5 = Ca*k_on_1*Prime_4 - (5*b^4*k_off_1P_plus*f^5)*Prime_5;
%     
%     dSuperprime_1 = 5*Ca*k_on_1*Superprime + 2*b*k_off_1*Superprime_2 - (k_off_1 + 4*Ca*k_on_1 + (M_plus + P_plus)*f)*Superprime_1;
%     dSuperprime_2 = 4*Ca*k_on_1*Superprime + 3*b^2*k_off_1*Superprime_3 - (2*b*k_off_1 + 3*Ca*k_on_1 + (M_plus + P_plus)*f^2)*Superprime_2;
%     dSuperprime_3 = 3*Ca*k_on_1*Superprime + 4*b^3*k_off_1*Superprime_4 - (3*b^2*k_off_1 + 2*Ca*k_on_1 + (M_plus + P_plus)*f^3)*Superprime_3;
%     dSuperprime_4 = 2*Ca*k_on_1*Superprime + 5*b^4*k_off_1*Superprime_5 - (4*b^3*k_off_1 + Ca*k_on_1 + (M_plus + P_plus)*f^4)*Superprime_4;
%     dSuperprime_5 = Ca*k_on_1*Superprime - (5*b^4*k_off_1 + (M_plus + P_plus)*f^5)*Superprime_5;
%     
%     
%     
%     dFused = M_plus*(Mature + f*Mature_1 + f^2*Mature_2 + f^3*Mature_3 + f^4*Mature_4 + f^5*Mature_5) + P_plus*(Prime + f*Prime_1 + f^2*Prime_2 + f^3*Prime_3 + f^4*Prime_4 + f^5*Prime_5) + (M_plus + P_plus)*(Superprime + f*Superprime_1 + f^2*Superprime_2 + f^3*Superprime_3 + f^4*Superprime_4 + f^5*Superprime_5);
%     
%     dCa = (-1/T_Ca_decay)*Ca + Ca_spike*ismember(t,stimulus_times);
%     
%     partials = [dImmature, dMature, dPrime, dSuperprime, dMature_1, dMature_3, dMature_4, dMature_5, dMature_5, dPrime_1, dFused, dCa];
% 
% end
% 