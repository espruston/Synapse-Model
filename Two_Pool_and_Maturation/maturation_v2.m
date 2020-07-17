%maturation_v2 assumes a maturation model where the maturation can be
%accelerated by calcium from a pulse

stimulus_times = linspace(20,200,10);
k_rep = .002;
k_mat = 8e8;
k_unmat = 120;
k_on_1 = 1.4e8;
k_off_1 = 4e3;
k_on_2 = 4e7;
k_off_2 = 2e3;
f = 28;
s = 510;
k_fuse_basal = 3.5e-4;
Ca_rest = 5e-8;
Ca_spike = 2.5e-5;
T_Ca_decay = 0.04;

m = 2; %second sensor Ca cooperativity
n = 5; %syt1 Ca cooperativity
b = 0.5; %sensor intercooperativity
n_sites = 1;

l = (m+1)*(n+1);

function partials = maturation_v2(t,y)

    state = y(1:l); %1:L
    immature = y(l+1);
    fused = y(l+2);
    Ca = y(l+3);
  
    
    
    dImmature_0 = ;
    dImmature_1 = ;
    dImmature_2 = ;
    dMature_0 = ;
    dMature_1 = ;
    dMature_2 = ;
    dMature_3 = ;
    dMature_4 = ;
    dMature_5 = ;
    
    
    
    dFused
    
    dCa = (-1/T_Ca_decay)*Ca + Ca_spike*ismember(t,stimulus_times);
    
    partials = [dFused, dCa];

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