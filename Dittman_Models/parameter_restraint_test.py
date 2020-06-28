import numpy as np
import pandas as pd
from math import e

if __name__ == "__main__":

    data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', nrows = 20)
    data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 22, nrows = 20)
    data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 44, nrows = 20)
    data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 66, nrows = 20)

    delta_D = 1
    delta_F = 1
    T_F = 100
    F_1 = 1-.103
    K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*0 - F_1) - 1)
    k_0 = 1
    #r = np.asarray([1,10,20,50])
    r = np.asarray([1, 10,20])

    T_D_valid = []
    K_D_valid = []
    k_max_valid = []

    n_pulses = 20

    ss_stdev_1hz = data_1hz.to_numpy()[n_pulses-1,1]*5
    ss_stdev_10hz = data_10hz.to_numpy()[n_pulses-1,1]*5
    ss_stdev_20hz = data_20hz.to_numpy()[n_pulses-1,1]*5
    ss_stdev_50hz = data_50hz.to_numpy()[n_pulses-1,1]*5

    ss_1hz = data_1hz.to_numpy()[n_pulses-1, 0]
    ss_10hz = data_10hz.to_numpy()[n_pulses-1, 0]
    ss_20hz = data_20hz.to_numpy()[n_pulses-1, 0]
    ss_50hz = data_50hz.to_numpy()[n_pulses-1, 0]

    ss_range_1hz = [ss_1hz - ss_stdev_1hz, ss_1hz + ss_stdev_1hz]
    ss_range_10hz = [ss_10hz - ss_stdev_10hz, ss_10hz + ss_stdev_10hz]
    ss_range_20hz = [ss_20hz - ss_stdev_20hz, ss_20hz + ss_stdev_20hz]
    ss_range_50hz = [ss_50hz - ss_stdev_50hz, ss_50hz + ss_stdev_50hz]


    #single hz test
    # for T_D in np.linspace(1, 1000, 1000):
    #     for K_D in np.linspace(0, 1000, 1001):
    #         for k_max in np.linspace(k_0, 99+k_0, 100):
    #             if 1/(r*T_F) > 709:
    #                 CaX_F_ss = 0
    #             else:
    #                 CaX_F_ss = delta_F*((e**(1/(r*T_F)) - 1)**(-1))
    #             if 1/(r*T_D) > 709:
    #                 CaX_D_ss = 0
    #             else:
    #                 CaX_D_ss = delta_D*((1 - e**(-1/(r*T_D)))**(-1))
    #
    #             #ss values used after iteration
    #             if CaX_F_ss == 0:
    #                 F_ss = F_1
    #             else:
    #                 F_ss = F_1 + (1 - F_1)/(1+(K_F/CaX_F_ss)) #eq 11
    #
    #             if CaX_D_ss == 0:
    #                 xi_ss = 1
    #             else:
    #                 xi_ss = ((K_D/CaX_D_ss + 1)/((K_D/CaX_D_ss) + e**(-1/(r*T_D))))**(-1*(k_max-k_0)*T_D) #modified eq 16
    #             D_ss = (1 - e**(-1*k_0/r)*xi_ss)/(1 - (1 - F_ss)*e**(-1*k_0/r)*xi_ss) #eq 20
    #             EPSC_norm_ss = D_ss*(F_ss/F_1) #eq 21
    #
    #             if ss_range_20hz[0] <= EPSC_norm_ss <= ss_range_20hz[1]: #all steady state values are within 5 SEM
    #                 T_D_valid.append(T_D)
    #                 K_D_valid.append(K_D)
    #                 k_max_valid.append(k_max)

    #multiple hz test
    for T_D in np.linspace(1, 1000, 1000):
        for K_D in np.linspace(0, 1000, 1001):
            for k_max in np.linspace(k_0, 99+k_0, 100):
                CaX_F_ss = np.zeros(len(r))
                CaX_F_ss[np.where(1/(r*T_F) < 709)] = delta_F*((e**(1/(r*T_F)) - 1)**(-1))
                CaX_D_ss = np.zeros(len(r))
                CaX_D_ss[np.where(1/(r*T_D) < 709)] = delta_D*((1 - e**(-1/(r*T_D)))**(-1))

                #ss values used after iteration
                F_ss = np.zeros(len(r))
                F_ss[np.where(CaX_F_ss == 0)] = F_1
                F_ss[np.where(CaX_F_ss != 0)] = F_1 + (1 - F_1)/(1+(K_F/CaX_F_ss)) #eq 11

                xi_ss = np.zeros(len(r))
                xi_ss[np.where(CaX_D_ss == 0)] = 1
                xi_ss[np.where(CaX_D_ss != 0)] = ((K_D/CaX_D_ss + 1)/((K_D/CaX_D_ss) + e**(-1/(r*T_D))))**(-1*(k_max-k_0)*T_D) #modified eq 16

                D_ss = (1 - e**(-1*k_0/r)*xi_ss)/(1 - (1 - F_ss)*e**(-1*k_0/r)*xi_ss) #eq 20
                EPSC_norm_ss = D_ss*(F_ss/F_1) #eq 21

                if (ss_range_1hz[0] <= EPSC_norm_ss[0] <= ss_range_1hz[1]) and (ss_range_10hz[0] <= EPSC_norm_ss[1] <= ss_range_10hz[1]) and (ss_range_20hz[0] <= EPSC_norm_ss[2] <= ss_range_20hz[1]):
                    T_D_valid.append(T_D)
                    K_D_valid.append(K_D)
                    k_max_valid.append(k_max)

    #all hz test
    # for T_D in np.linspace(1, 1000, 1000):
    #     for K_D in np.linspace(0, 1000, 1001):
    #         for k_max in np.linspace(k_0, 99+k_0, 100):
    #             CaX_F_ss = np.zeros(4)
    #             CaX_F_ss[np.where(1/(r*T_F) < 709)] = delta_F*((e**(1/(r*T_F)) - 1)**(-1))
    #             CaX_D_ss = np.zeros(4)
    #             CaX_D_ss[np.where(1/(r*T_D) < 709)] = delta_D*((1 - e**(-1/(r*T_D)))**(-1))
    #
    #             #ss values used after iteration
    #             F_ss = np.zeros(4)
    #             F_ss[np.where(CaX_F_ss == 0)] = F_1
    #             F_ss[np.where(CaX_F_ss != 0)] = F_1 + (1 - F_1)/(1+(K_F/CaX_F_ss)) #eq 11
    #
    #             xi_ss = np.zeros(4)
    #             xi_ss[np.where(CaX_D_ss == 0)] = 1
    #             xi_ss[np.where(CaX_D_ss != 0)] = ((K_D/CaX_D_ss + 1)/((K_D/CaX_D_ss) + e**(-1/(r*T_D))))**(-1*(k_max-k_0)*T_D) #modified eq 16
    #
    #             D_ss = (1 - e**(-1*k_0/r)*xi_ss)/(1 - (1 - F_ss)*e**(-1*k_0/r)*xi_ss) #eq 20
    #             EPSC_norm_ss = D_ss*(F_ss/F_1) #eq 21

                # if (ss_range_1hz[0] <= EPSC_norm_ss[0] <= ss_range_1hz[1]) and (ss_range_10hz[0] <= EPSC_norm_ss[1] <= ss_range_10hz[1]) and (ss_range_20hz[0] <= EPSC_norm_ss[2] <= ss_range_20hz[1]) and (ss_range_50hz[0] <= EPSC_norm_ss[3] <= ss_range_50hz[1]): #all steady state values are within 5 SEM
                #     T_D_valid.append(T_D)
                #     K_D_valid.append(K_D)
                #     k_max_valid.append(k_max)
