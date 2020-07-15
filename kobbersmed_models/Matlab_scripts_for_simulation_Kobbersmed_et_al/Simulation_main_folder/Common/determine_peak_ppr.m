function [peak1, peak2, peaksum, ppr, peak2_extr, ppr_extr, ttp2080, peaksum_extr] = determine_peak_ppr(time_vector, EPSC, stim_freq)
    
plot_on_off = 0;
expo_fit_on_off = 1;

ISI = 1000/stim_freq;

t_ISI_ind = find(time_vector >= (0.5+ISI)*1e-3, 1, 'first');
% t_15ms_ind = find(time_vector >= 15e-3, 1, 'first');

[peak1, peak1_ind] = min(EPSC(1:t_ISI_ind)); %peak1 and index

[peak2_0, peak2_ind_0] = min(EPSC(t_ISI_ind:end)); %peak2 measured from baseline

peak2_ind = peak2_ind_0 + t_ISI_ind - 1; %peak2 index
    
if peak1 ~= 0
    peak20_ind = find(EPSC<=(0.2*peak1),1,'first'); %Location of 20% of peak
    peak80_ind = find(EPSC<=(0.8*peak1),1,'first'); %location of 80% of peak
    [cleft, cleft_ind_0] = max(EPSC(peak1_ind:peak2_ind));
    cleft_ind = cleft_ind_0 + peak1_ind - 1;
    peak2 = peak2_0-cleft;
    ppr = peak2/peak1;
    peaksum = peak1+peak2;
    ttp2080 = time_vector(peak80_ind) - time_vector(peak20_ind);
    
    %%%%Extrapolation: Peak 2 is determined by estimating
    %%%%exponential decay of the first response.
    
    fit_ind = peak1_ind:cleft_ind;
    
    fit_ind2 = (peak1_ind+round(0.5*length(fit_ind))):cleft_ind;
    
    
    

    if expo_fit_on_off == 1    
        if length(fit_ind) > 1
            expo_fit = fit(time_vector(fit_ind), EPSC(fit_ind), 'exp1');
            expo_fit_alt = fit(time_vector(fit_ind2), EPSC(fit_ind2), 'exp1');
            peak2BaseExtr = expo_fit(time_vector(peak2_ind));
            peak2BaseExtr_alt = expo_fit_alt(time_vector(peak2_ind));
            peak2_extr = peak2_0 - peak2BaseExtr;
            ppr_extr = peak2_extr/peak1;
            peaksum_extr = peak1 + peak2_extr;
        else
            [peaksum_extr, peak2_extr, ppr_extr]  = deal(NaN);
        end
    else
            [peaksum_extr, peak2_extr, ppr_extr]  = deal(NaN);
    end

elseif peak1 == 0 || issorted(EPSC(1e4:2e4))
    [peak2, peaksum, ppr, peak2_extr, ppr_extr, ttp2080, peaksum_extr] = deal(NaN);
end


%%%Illustration of the extrapolation

if plot_on_off ==1
    figure
    hold on
    plot(time_vector*1e3, EPSC*1e9, 'b', 'LineWidth', 2)
    plot(time_vector*1e3, expo_fit(time_vector)*1e9, 'r', 'LineWidth', 1)
    plot([time_vector(peak2_ind) time_vector(peak2_ind)]*1e3, [peak2BaseExtr EPSC(peak2_ind)]*1e9, 'k--')
    xlim([0 20])
    ylim([min(EPSC)*1e9 0])
    xlabel('Time [ms]')
    ylabel('Current [nA]')
    set(gca, 'xtick', [0 10 20], 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out')

    figure
    hold on
    plot(time_vector*1e3, EPSC*1e9, 'b', 'LineWidth', 2)
    plot(time_vector*1e3, expo_fit_alt(time_vector)*1e9, 'r', 'LineWidth', 1)
    plot([time_vector(peak2_ind) time_vector(peak2_ind)]*1e3, [peak2BaseExtr_alt EPSC(peak2_ind)]*1e9, 'k--')
    xlim([0 20])
    ylim([min(EPSC)*1e9 0])
    xlabel('Time [ms]')
    ylabel('Current [nA]')
    set(gca, 'xtick', [0 10 20], 'LineWidth', 1, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out')
end


