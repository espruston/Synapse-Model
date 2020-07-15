function explore_Caprim_results

CaExtracellular = [0.75 1.5 3 6 10];
Ca_prim_type = 3;%0;
% prim_kMs = [60e-9 80e-9 100e-9 120e-9];%0;
prim_kMs = [65e-9 70e-9 75e-9 100e-9 120e-9 150e-9 200e-9];
prim_rate_consts = [100 150 180 210 250 300 500];%0;
unprim_rate_consts = [50 100 200 300 400 500];%0;

Q_max = 12;

model_type = 41;

size_results = [length(prim_kMs) length(prim_rate_consts) length(unprim_rate_consts)];
length_results = prod(size_results);
num_calc = length(CaExtracellular);

pprs_prim = zeros(length_results, num_calc);
peak1_prim = zeros(length_results, num_calc);

for k = 1:length(prim_kMs)
    for l = 1:length(prim_rate_consts)
        for m = 1:length(unprim_rate_consts)

            prim_kM = prim_kMs(k);
            prim_rate_const = prim_rate_consts(l); 
            unprim_rate_const = unprim_rate_consts(m);
            
            par_free = [Q_max Ca_prim_type prim_kM prim_rate_const unprim_rate_const];

            [par_init, savefilename] = parameter_choices(par_free, model_type)

            filename = [savefilename '_DET_fewresults.mat']
            
            load(filename)
            
            ind = sub2ind(size_results, k, l, m);
            
            pprs_prim(ind,:) = pprs_extr;
            peak1_prim(ind,:) = peak1s;
        end
    end
end


load('DataSummary_extr.mat')

ppr_exp = repmat(ppr_extr_mean_all_cells_exp, length_results,1);
peak1_exp = repmat(peak1_mean_all_cells_exp, length_results,1);

cost = sum( ( ((ppr_exp - pprs_prim).^2) ./ ppr_exp ));% + ( ((peak1_exp - peak1_prim).^2) ./ abs(peak1_exp)),2);

[costmin, costmin_ind] = min(cost);

[km_ind, prim_ind, unprim_ind] =  ind2sub(size_results, costmin_ind);

costmin_ind = 99;


best_kM = prim_kMs(km_ind);
best_prim_rate = prim_rate_consts(prim_ind);
best_unprim_rate = unprim_rate_consts(unprim_ind);

best_ppr = pprs_prim(costmin_ind,:);
best_peak1 = peak1_prim(costmin_ind,:);


figure
hold on

yyaxis left
errorbar(CaExtracellular, abs(peak1_mean_all_cells_exp), peak1_std_of_means_exp, 'o', 'Color', 'k', 'MarkerEdgeColor', 'k')
% plot(CaExtracellular, abs(peak1s_notcorr), '-o', 'LineWidth', 2, 'Color', 'b', 'MarkerFacecolor', 'b')
plot(CaExtracellular, abs(best_peak1), 'o', 'Color', 'r')

xlabel('Extracellular calcium [mM]')
xlim([0 11])
ylabel('Absolute peak [nA]')
ylim([0 120])
set(gca, 'xtick', [0 5 10], 'ytick', [0 50 100])

yyaxis right

ax = gca;
ax.ColorOrderIndex = 1;


errorbar(CaExtracellular, ppr_extr_mean_all_cells_exp, ppr_std_of_means_exp, '^', 'Color', 'k')
% plot(CaExtracellular, pprs_extr_notcorr, '--o', 'LineWidth', 2, 'Color', 'b', 'MarkerFacecolor', 'b')
plot(CaExtracellular, best_ppr, '^', 'Color', 'r')
plot([0 15], [1 1])
ylabel('ppr')
ylim([0 2.3])
set(gca, 'ytick', [0 1 2], 'TickDir', 'out')

legend('Peak (experiment)', 'Peak (simulation)', 'ppr (experiment)', 'ppr (simulation)')


saveas(gcf, './Sim_data/new_figures/Caprim_exploration.eps', 'epsc')





%%%%%If all states are saved: 



if (exist('states_all_cell)')>0
    
            par_free = [Q_max Ca_prim_type best_kM prim_rate_const unprim_rate_const];

            [par_init, savefilename] = parameter_choices(par_free, model_type)

            filename = [savefilename '_DET_fewresults.mat']
            
            load(filename)
            

            
            figure
            hold on


end