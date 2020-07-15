function [fast_states, slow_states, total_fast_states, total_slow_states] = determine_total_slowfast_DET(time_vector, states_all, par)


n_max = par(13);
m_max = par(14);

fast_states = zeros(length(time_vector), n_max + 1);
slow_states = zeros(length(time_vector), m_max + 1);
total_fast_states = zeros(length(time_vector),1);
total_slow_states = zeros(length(time_vector),1);

for k = 1:(n_max+1)
    fast_states(:,k) = sum(states_all(:,((k-1)*(m_max+1)+1):(k*(m_max+1))),2); %Fast states changing over time
    total_fast_states = total_fast_states + (k-1)*fast_states(:,k); %Total number of ions bound
end

for k = 1:(m_max+1)
    slow_states(:,k) = sum(states_all(:,k:(m_max+1):(n_max*(m_max+1)+k)),2); %Slow states changing over time
    total_slow_states = total_slow_states + (k-1)*slow_states(:,k); %Total number of ions bound
end
