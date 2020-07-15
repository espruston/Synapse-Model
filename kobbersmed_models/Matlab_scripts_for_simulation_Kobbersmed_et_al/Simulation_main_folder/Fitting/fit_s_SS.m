function [opt_s, cost_total] = fit_s_SS(Q_max)


k_4 = (1.4e8)/2;
k_min4 = 4000/60;



filecont = sprintf(     ['-------------------' '\n' ...
                         '-------------------' '\n' ...
                         'Now Optimizing to the following values:' '\n' ...
                         'Q_max: ' num2str(Q_max) '\n' ...
                         'k_4: ' num2str(k_4) '\n' ...
                         'k_min4: ' num2str(k_min4) '\n' ...
                         '------------------- \n', ...
                         '------------------- \n' '\n' ...
                        ]);

fullfile = ['./Sim_data/Log_files/Log_Optimization_s_SS'];

newParameterfile = fopen(fullfile,'a');
fprintf(newParameterfile,'%s\n',filecont);
fclose(newParameterfile);

fclose('all');



s_init = 100;

cost_func = @(s_val)cost_value_SS_Q_max(Q_max, k_4, k_min4, s_val);

[opt_s, cost_total] = fminsearch(cost_func, s_init);



