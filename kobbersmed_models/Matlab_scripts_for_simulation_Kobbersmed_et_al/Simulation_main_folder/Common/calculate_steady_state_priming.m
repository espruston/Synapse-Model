function steady_priming = calculate_steady_state_priming(par, Calcium)


a =par(32);
b =par(33);

k_on =par(6);
k_off =par(7);

u = par(43);
m_max = par(14);
num_ves = par(22);

steady_priming_nonorm = zeros(length(Calcium),m_max+2);
steady_priming = zeros(length(Calcium),m_max+2);

for l = 1:length(Calcium)
    Ca = Calcium(l);
steady_priming_nonorm(l,4) = 100;
steady_priming_nonorm(l,3) = ((2*k_off + b*u^(-2))./(k_on*Ca) )*steady_priming_nonorm(l,4);
steady_priming_nonorm(l,2) = ((k_off + k_on*Ca + b*u^(-1))*steady_priming_nonorm(l,3) - 2*k_off*steady_priming_nonorm(l,4))./(2*k_on*Ca);
steady_priming_nonorm(l,1) = sum((b*steady_priming_nonorm(l,2:end).*u.^(-(0:m_max))),2)/a;

steady_priming(l,:) = steady_priming_nonorm(l,:)/sum(steady_priming_nonorm(l,:));
end


priming1 = 1-steady_priming(:,1);