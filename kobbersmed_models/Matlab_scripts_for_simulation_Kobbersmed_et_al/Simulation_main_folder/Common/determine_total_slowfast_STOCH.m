function [fast_states, slow_states, total_fast_states, total_slow_states] = determine_total_slowfast_STOCH(ves_states, par)

n_max = par(13);
m_max = par(14);

ves_states_mod = ves_states; %Remove all fusion/empty
inds_mor100 = (ves_states>= 100);
ves_states_mod(inds_mor100) = ves_states_mod(inds_mor100)-100;

FS = ves_states_mod - floor(ves_states_mod/10);
SS = floor(ves_states_mod/10);

total_fast_states = sum(FS);
total_slow_states = sum(SS);

states_size = size(ves_states);
length_states = states_size(2);

fast_states = zeros(length_states, n_max + 1);
slow_states = zeros(length_states, m_max + 1);


for k = 1:(n_max+1)
    fast_states(:,k) = sum(FS == k-1 ); %Fast states changing over time
end

for k = 1:(m_max+1)
    slow_states(:,k) = sum(SS == k-1);
end

% fast_total_test = zeros(6246,1);
% for k = 1:(n_max+1)
%     fast_total_test = fast_total_test + (k-1)*(fast_states(:,k));
% end
% 


