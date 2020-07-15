function [activation_act0_norep, delay_act0_norep, inactivation_act0_norep] = determine_activation_from_simdata(time_vector, states_all, plot_on_off, n_max, m_max)


no_of_Ca = 1;

activation_act0_norep = zeros(length(time_vector), no_of_Ca);
delay_act0_norep = zeros(length(time_vector), no_of_Ca);
inactivation_act0_norep = zeros(length(time_vector), no_of_Ca);

for k = 1:no_of_Ca; %1:5
%     states_all = states_all_cell_act0_norep{k};
    total_active = sum([states_all(:,1:6)  states_all(:,19)], 2);
    total_delay = sum([states_all(:,7:12)  states_all(:,20)], 2);
    total_inactive = sum([states_all(:,13:18)  states_all(:,21)], 2);
    activation_act0_norep(:,k) = total_active./(total_inactive+ total_delay + total_active);
    delay_act0_norep(:,k) = total_delay./(total_inactive+ total_delay + total_active);
    inactivation_act0_norep(:,k) = total_inactive./(total_inactive+ total_delay + total_active);
end

% par = model_parameters_det(par_int, 14, 250, 0, act_model_type);

% k_M_act = par(25);

% Ca_activation_residual = (Ca_residuals./(Ca_residuals + k_M_act));
%   [~, ~, dist_8ms, Ca_8ms] = read_and_reshape_residual(file_to_read);
% 
%     min_activation = (Ca_8ms(100)./(Ca_8ms(100) + k_M_act));
%     max_activation = (Ca_8ms(100)./(Ca_8ms(100) + k_M_act));


if plot_on_off
    figure
    plot(time_vector*1e3, activation_act0_norep, 'LineWidth', 2)%, 'LineStyle', '--')
    hold on
    plot(time_vector*1e3, delay_act0_norep, 'LineWidth', 2)%, 'LineStyle', '--')
    plot(time_vector*1e3, inactivation_act0_norep, 'LineWidth', 2)%, 'LineStyle', '--')

            ax = gca;
            ax.ColorOrderIndex = 1;
    % plot([time_vector(1)*1e3 time_vector(20501)*1e3], [Ca_activation_residual; Ca_activation_residual], '--')
    xlabel('Time [ms]')
    ylabel('Fraction of sites')
    ylim([0 1])
    legend('A state', 'D state', 'I state')
end




