function mEPSC_fit()

load('mEPSC_new.mat')
load('mEPSC_aligned.mat')


time_vec_gen = 0:1e-6:34*1e-3;


A = -2.254687965059888e-06;
B = 0.5; %-37.227291574661400;
t_0 = 0.002498798993959;
tau_rf = 125.8688851183726;
tau_df = 8.440530570966260e-04;
tau_ds = 0.002946178296457;

pars = [A B t_0 tau_rf tau_df tau_ds];

numpars = length(pars);

cost_func1 =  @(par_int)mEPSC_cost1(par_int);
% cost_func2 =  @(par_int)mEPSC_cost2(par_int);
cost_func3 =  @(par_int)mEPSC_cost3(par_int);

fit_options = optimset('MaxFunEvals', 1e5*numpars, 'MaxIter', 1e5*numpars);

[best_par1, best_cost1_first_new] = fminsearch(cost_func1, pars, fit_options)
pars = best_par1;
[best_par1, best_cost1_second_new] = fminsearch(cost_func1, pars, fit_options)

A1 = best_par1(1);
B1 = best_par1(2);
t_01 = best_par1(3);
tau_rf1 = best_par1(4);
tau_df1 = best_par1(5);
tau_ds1 = best_par1(6);

mEPSC_gen_func1 = @(t)(t>=t_01).*(A1*(1-exp(-(t-t_01)/tau_rf1)).*(B1*exp(-(t-t_01)/tau_df1) + (1-B1).*exp(-(t-t_01)/tau_ds1)));

mEPSC_gen1 = mEPSC_gen_func1(time_vec_gen);


% 
% pars = best_par1;
% 
% [best_par2, best_cost2_first_new] = fminsearch(cost_func2, pars, fit_options)
% pars = best_par2;
% [best_par2, best_cost2_second_new] = fminsearch(cost_func2, pars, fit_options)
% 
% A2 = best_par2(1);
% B2 = best_par2(2);
% t_02 = best_par2(3);
% tau_rf2 = best_par2(4);
% tau_df2 = best_par2(5);
% tau_ds2 = best_par2(6);
% 
% mEPSC_gen_func2 = @(t)(t>=t_02).*(A2*(1-exp(-(t-t_02)/tau_rf2)).*(B2*exp(-(t-t_02)/tau_df2) + (1-B2).*exp(-(t-t_02)/tau_ds2)));
% 
% mEPSC_gen2 = mEPSC_gen_func2(time_vec_gen);


% A3 = -2.254687965059888e-06;
% B3 = -37.227291574661400;
% t_03 = 0.002498798993959;
% tau_rf3 = 125.8688851183726;
% tau_df3 = 8.440530570966260e-04;
% tau_ds3 = 0.002946178296457;

pars = best_par1;
pars(2) = 0.5

[best_par3, best_cost3_first_new] = fminsearch(cost_func3, pars, fit_options)
pars = best_par3;
[best_par3, best_cost3_second_new] = fminsearch(cost_func3, pars, fit_options)

A3 = best_par3(1);
B3 = best_par3(2);
t_03 = best_par3(3);
tau_rf3 = best_par3(4);
tau_df3 = best_par3(5);
tau_ds3 = best_par3(6);



mEPSC_gen_func3 = @(t)(t>=t_03).*(A3*(1-exp(-(t-t_03)/tau_rf3)).*(B3*exp(-(t-t_03)/tau_df3) + (1-B3).*exp(-(t-t_03)/tau_ds3)));

mEPSC_gen3 = mEPSC_gen_func3(time_vec_gen);



A = [best_par1(1) best_par3(1)]
B = [best_par1(2) best_par3(2)]
t_0 = [best_par1(3) best_par3(3)]
tau_r = [best_par1(4) best_par3(4)]
tau_df = [best_par1(5) best_par3(5)]
tau_ds = [best_par1(6) best_par3(6)]

system('rm ./Sim_data/mEPSC_fit_table.xlsx')

% VarNames = {'A'; 'Fit_2'; 'Fit_3'};
RowNames = {'Fit_1'; 'Fit_3'};
VarNames = {'A', 'B', 't_0', 'tau_r', 'tau_df', 'tau_ds'};
%Table: K_rep, s, f, k_d_first, k_3, k_d_second, k_4, pprs, cost
mEPSC_fit_table = table(A', B', t_0', tau_r', tau_df', tau_ds', 'RowNames', RowNames, 'VariableNames', VarNames);


writetable(mEPSC_fit_table, './Sim_data/mEPSC_fit_table.xlsx', 'WriteRowNames',true)






figure%('visible', 'off')
hold on
plot(mEPSC_time*1e3, mEPSC_current*1e9, 'color', 'k', 'LineWidth', 2)
plot(mEPSC_time*1e3, mEPSC_ave_final*1e9, 'color', 'b', 'LineWidth', 2)
% plot(time_vec_gen*1e3, mEPSC_gen2*1e9, 'LineWidth', 2, 'color', 'b')
plot(time_vec_gen*1e3, mEPSC_gen3*1e9, 'LineWidth', 2, 'color', 'r')
plot(mEPSC_time_smooth*1e3, mEPSC_current_smooth*1e9, '--', 'color', 'k', 'LineWidth', 2)
legend('Old experimental mini', 'New experimental mini', 'Best fit')
xlabel('Time [ms]')
ylabel('Current [nA]')
set(gca, 'TickDir', 'out')
saveas(gcf, './Ad_hoc_figures/Meeting_190304/newmini_fit.eps', 'epsc')

figure('visible', 'off')
hold on
plot(mEPSC_time*1e3, mEPSC_current*1e9, 'color', 'k', 'LineWidth', 2)
plot(time_vec_gen*1e3, mEPSC_gen1*1e9, 'LineWidth', 2, 'color', 'b')
legend('Old experimental mini', 'Best fit')
xlabel('Time [ms]')
ylabel('Current [nA]')
set(gca, 'TickDir', 'out')
saveas(gcf, './Ad_hoc_figures/Meeting_190304/oldmini_fit.eps', 'epsc')

