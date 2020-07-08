import numpy as np
import pandas as pd
from TwoPoolAndMat import two_pool
from matplotlib.widgets import Slider
from matplotlib import pyplot as plt
from sklearn import metrics
from math import e

if __name__ == "__main__":

    data_WT = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Raw data climbing fiber bursts DW.xlsx', usecols = 'A,B,C', skiprows = 93, nrows = 11)

    ITIs = data_WT.to_numpy()[:,0]
    print(ITIs)
    EPSCs2_WT = data_WT.to_numpy()[:,1]
    stdev_WT = data_WT.to_numpy()[:,2]

    r = 50

    m = .889
    size_fast = m
    size_slow = 1-size_fast

    p_fast = 0.56
    p_slow = 1-.198

    T_fast = .2 #in seconds
    T_slow = 7

    F = 1
    K_F = 1.5
    delta_F = 1
    T_F = 0.0001
    BG_F = 1e-5
    sat_F = K_F
    facil = BG_F

    n_pulses = 20

    first_burst = two_pool(50, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F)

    fastpool_initial = size_fast
    slowpool_initial = size_slow
    EPSCs2_reg = []

    for t in ITIs:
        state = first_burst.final_state
        delta_t = t - 20 #1 ISI has already been simulated for the final state
        r = 1000/delta_t

        fastrecover = (fastpool_initial - state[0])*(1-e**(-1/(r*T_fast)))
        state[0] += fastrecover
        slowrecover = (slowpool_initial - state[1])*(1-e**(-1/(r*T_slow)))
        state[1] += slowrecover

        EPSCs2_reg.append((state[0]*p_fast*F + state[1]**p_slow*F)/first_burst.norm_val)

    #
    # best = [-1e30]
    # search_p_fast = np.linspace(.01,1,10)
    # search_T_fast = np.linspace(.1,1,10)
    # search_T_slow = np.linspace(1,10,10)
    # search_size_fast = np.linspace(0,1,10)
    #
    # n_consider = 20

    # for p_fast in search_p_fast:
    #     for T_fast in search_T_fast:
    #         for T_slow in search_T_slow:
    #             for size_fast in search_size_fast:
    #
    #                 size_slow = 1 - size_fast
    #
    #                 EPSC_1 = two_pool(1, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
    #                 EPSC_10 = two_pool(10, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
    #                 EPSC_20 = two_pool(20, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
    #                 EPSC_50 = two_pool(50, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
    #
    #                 r_squared_1hz = metrics.r2_score(EPSC_1hz[0:n_consider], EPSC_1[0:n_consider])
    #                 r_squared_10hz = metrics.r2_score(EPSC_10hz[0:n_consider], EPSC_10[0:n_consider])
    #                 r_squared_20hz = metrics.r2_score(EPSC_20hz[0:n_consider], EPSC_20[0:n_consider])
    #                 r_squared_50hz = metrics.r2_score(EPSC_50hz[0:n_consider], EPSC_50[0:n_consider])
    #
    #                 avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4
    #
    #                 if avg_r_squared > best[-1]:
    #                     best = [EPSC_1, EPSC_10, EPSC_20, EPSC_50, [p_fast, T_fast, T_slow, size_fast], [r_squared_1hz, r_squared_10hz, r_squared_20hz, r_squared_50hz], avg_r_squared]

    plt.errorbar(ITIs, EPSCs2_WT, yerr = stdev_WT, label = "data")
    plt.plot(ITIs, EPSCs2_reg, label = "regression")
    print(metrics.r2_score(EPSCs2_WT, EPSCs2_reg))

    # print('Best R squared:', best[-1])
    # print(best[-2])
    plt.legend()
    plt.show()
