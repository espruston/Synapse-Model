function [iChannel] = calculate_ichan(Q_max, K_m, CaExtracellular,  hill)

iChannel = (Q_max*(CaExtracellular.^hill))./((CaExtracellular.^hill)+(K_m^hill));

