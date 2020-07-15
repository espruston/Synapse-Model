function [ves_avail] = determine_ves_avail(time_vector, states_all, n_max, m_max)



num_states = (n_max + 1) * (m_max + 1);

num_ves_avail = sum(states_all(:,1:num_states), 2);

ves_avail = num_ves_avail / sum(states_all(1,:));


