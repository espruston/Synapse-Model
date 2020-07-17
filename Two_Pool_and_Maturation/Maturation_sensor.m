%maturation_v2 assumes a maturation model where the maturation can be
%accelerated by calcium from a pulse

k_refill = 2; %vesicles ms-1

k_mature = 8e8; %
k_unmature = 120; %

k_on_1 = 1.4e5; %ca binding rate, syt1, M-1ms-1
k_off_1 = 4; %ca unbinding rate, syt1, ms-1 
%k_on_2 = 4e4; %ca binding rate, syt3, M-1ms-1, estimate
%k_off_2 = .5; %ca unbinding rate syt3, ms-1, estimate

b = 0.5; %cooperativity of Ca binding syt1
%c = 1; %cooperativity of syt3 for staying matured

f = 28;
k_fuse_basal = 3.5e-4;
M_plus = 3.5e-4;
P_plus = 3.5e-4;

Ca_rest = 5e-8;
Ca_spike = 2.5e-5;
Ca_residual = 250e-9;
T_Ca_decay = 0.04;

max_time = 500;
stimulus_times = [0,20];
delta_t = 1e-3;

FWHM = .34; %Local calcium full width half maximum ms
sigma = FWHM/2.35; %variance
mu = 2*FWHM; %time at which Ca_spike is maximal (ms)

ts = linspace(0,max_time,max_time/delta_t);

Ca_sim = zeros(length(ts));
Ca_sim = Ca_sim + Ca_rest;

for t = 1:length(stimulus_times)
   
    spike_start_index = stimulus_times(t)/delta_t;
    spike_peak_index = (stimulus_times(t)+mu)/delta_t;
    
    Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end-spike_start_index) - mu)/sigma).^2/2);
    
    Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index)/T_Ca_decay);    

end

state_ini = [0,1,0,0]; %immature, mature, Prime, and Superprime states


for t = 1:length(ts)
   
    Ca = Ca_sim(t);
    
    
    
    
    
end

function partials = maturation_v2(t,y)

    Immature = y(1);
    Mature = y(2);
    Prime = y(3);
    Superprime = y(4);
    
    Mature_1 = y(5);
    Mature_2 = y(6);
    Mature_3 = y(7);
    Mature_4 = y(8);
    Mature_5 = y(9);
    
    Prime_1 = y(10);
    Prime_2 = y(11);
    Prime_3 = y(12);
    Prime_4 = y(13);
    Prime_5 = y(14);
    
    Superprime_1 = y(15);
    Superprime_2 = y(16);
    Superprime_3 = y(17);
    Superprime_4 = y(18);
    Superprime_5 = y(19);
    
    Fused = y(20);
    Ca = y(21);
  
    
    
    dImmature_ = k_refill*(1-Immature) - (Ca*k_mature + Ca*k_prime)*Immature + k_unmature*Mature + k_unprime*Prime;
    dMature = Ca*k_mature + k_unprime*Superprime + k_off_1*Mature_1 - (k_immature + Ca*k_prime + 5*Ca*k_on_1 + M_plus)*Mature;
    dPrime = Ca*k_prime*Immature + k_unprime*Superprime + k_off_1*Prime_1 - (k_unprime + Ca*k_mature + 5*Ca*k_on_1 + P_plus)*Prime;
    dSuperprime = Ca*k_prime*Mature + Ca*k_mature*Prime + k_off_1*Superprime_1 - (k_unprime + k_unmature + 5*Ca*k_on_1 + (M_plus + P_plus))*Superprime;
    
    dMature_1 = 5*Ca*k_on_1*Mature + 2*b*k_off_1*Mature_2 - (k_off_1 + 4*Ca*k_on_1 + M_plus*f)*Mature_1;
    dMature_2 = 4*Ca*k_on_1*Mature_1 + 3*b^2*k_off_1*Mature_3 - (2*b*k_off_1 + 3*Ca*k_on_1 + M_plus*f^2)*Mature_2;
    dMature_3 = 3*Ca*k_on_1*Mature_2 + 4*b^3*k_off_1*Mature_4 - (3*b^2*k_off_1 + 2*Ca*k_on_1 + M_plus*f^3)*Mature_3;
    dMature_4 = 2*Ca*k_on_1*Mature_3 + 5*b^4*k_off_1*Mature_5 - (4*b^3*k_off_1 + Ca*k_on_1 + M_plus*f^4)*Mature_4;
    dMature_5 = Ca*k_on_1*Mature_4 - (5*b^3*k_off_1 + M_plus*f^5)*Mature_5;
    
    dPrime_1 = 5*Ca*k_on_1*Prime + 2*b*k_off_1*Prime_2 - (k_off_1 + 4*Ca*k_on_1 + P_plus*f)*Prime_1;
    dPrime_2 = 4*Ca*k_on_1*Prime_1 + 3*b^2*k_off_1*Prime_3 - (2*b*k_off_1 + 3*Ca*k_on_1 + P_plus*f^2)*Prime_2;
    dPrime_3 = 3*Ca*k_on_1*Prime_2 + 4*b^3*k_off_1*Prime_4 - (3*b^2*k_off_1 + 2*Ca*k_on_1 + P_plus*f^3)*Prime_3;
    dPrime_4 = 2*Ca*k_on_1*Prime_3 + 5*b^4*k_off_1*Prime_5 - (4*b^3*k_off_1 + Ca*k_on_1 + P_plus*f^4)*Prime_4;
    dPrime_5 = Ca*k_on_1*Prime_4 - (5*b^4*k_off_1P_plus*f^5)*Prime_5;
    
    dSuperprime_1 = 5*Ca*k_on_1*Superprime + 2*b*k_off_1*Superprime_2 - (k_off_1 + 4*Ca*k_on_1 + (M_plus + P_plus)*f)*Superprime_1;
    dSuperprime_2 = 4*Ca*k_on_1*Superprime + 3*b^2*k_off_1*Superprime_3 - (2*b*k_off_1 + 3*Ca*k_on_1 + (M_plus + P_plus)*f^2)*Superprime_2;
    dSuperprime_3 = 3*Ca*k_on_1*Superprime + 4*b^3*k_off_1*Superprime_4 - (3*b^2*k_off_1 + 2*Ca*k_on_1 + (M_plus + P_plus)*f^3)*Superprime_3;
    dSuperprime_4 = 2*Ca*k_on_1*Superprime + 5*b^4*k_off_1*Superprime_5 - (4*b^3*k_off_1 + Ca*k_on_1 + (M_plus + P_plus)*f^4)*Superprime_4;
    dSuperprime_5 = Ca*k_on_1*Superprime - (5*b^4*k_off_1 + (M_plus + P_plus)*f^5)*Superprime_5;
    
    
    
    dFused = M_plus*(Mature + f*Mature_1 + f^2*Mature_2 + f^3*Mature_3 + f^4*Mature_4 + f^5*Mature_5) + P_plus*(Prime + f*Prime_1 + f^2*Prime_2 + f^3*Prime_3 + f^4*Prime_4 + f^5*Prime_5) + (M_plus + P_plus)*(Superprime + f*Superprime_1 + f^2*Superprime_2 + f^3*Superprime_3 + f^4*Superprime_4 + f^5*Superprime_5);
    
    dCa = (-1/T_Ca_decay)*Ca + Ca_spike*ismember(t,stimulus_times);
    
    partials = [dImmature, dMature, dPrime, dSuperprime, dMature_1, dMature_3, dMature_4, dMature_5, dMature_5, dPrime_1, dFused, dCa];

end


state_ini = zeros(1,l+3);
c = 0;
for i = 0:m+1
   for j = 0:n+1
       state_ini[c] = (((factorial(n)/factorial(n - j))*(Ca_rest.^j)*(k_on_1.^j))/(factorial(j)*(b.^(j*(j-1)/2))*(k_off_1.^j)))*(((factorial(m)/factorial(m - i))*(Ca_rest.^i)*(k_on_2.^i))/(factorial(i)*(b.^(i*(i-1)/2))*(k_off_2.^i)))
       c = c + 1;
   end
end
state_ini(l+1) = n_sites - sum(state_ini(0:l+1)); %number of immature vesicles
state_ini(l+2) = 0; %number of fused vesicles
state_ini(l+3) = Ca_rest;




max_time = stimulus_times(length(stimulus_times)) + stimulus_times(0)*2;
tspan = [0 max_time];
ts = linspace(0, max_time, max_time+1);

sol = ode45(@partials, tspan, state_ini);