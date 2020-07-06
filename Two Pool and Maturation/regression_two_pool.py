import numpy as np
import pandas as pd
import os
from TwoPoolAndMat import two_pool
from matplotlib.widgets import Slider
from matplotlib import pyplot as plt
from sklearn import metrics

if __name__ == "__main__":

    data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', nrows = 20)
    data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 22, nrows = 20)
    data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 44, nrows = 20)
    data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 66, nrows = 20)

    EPSC_1hz = data_1hz.to_numpy()[:,0]
    stdev_1hz = data_1hz.to_numpy()[:,1]
    EPSC_10hz = data_10hz.to_numpy()[:,0]
    stdev_10hz = data_10hz.to_numpy()[:,1]
    EPSC_20hz = data_20hz.to_numpy()[:,0]
    stdev_20hz = data_20hz.to_numpy()[:,1]
    EPSC_50hz = data_50hz.to_numpy()[:,0]
    stdev_50hz = data_50hz.to_numpy()[:,1]

    m = .68
    size_fast = m
    size_slow = 1-size_fast

    p_fast = 0.02
    p_slow = 1-.198

    T_fast = .11
    T_slow = 4.9

    F = 1
    K_F = 1.5
    delta_F = 1
    T_F = 0.0001
    BG_F = 1e-5
    sat_F = K_F
    facil = BG_F

    n_pulses = 20

    best = [-1e30]
    search_p_fast = np.linspace(.01,1,10)
    search_T_fast = np.linspace(.1,1,10)
    search_T_slow = np.linspace(1,10,10)
    search_size_fast = np.linspace(0,1,10)

    n_consider = 20

    for p_fast in search_p_fast:
        for T_fast in search_T_fast:
            for T_slow in search_T_slow:
                for size_fast in search_size_fast:

                    size_slow = 1 - size_fast

                    EPSC_1 = two_pool(1, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
                    EPSC_10 = two_pool(10, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
                    EPSC_20 = two_pool(20, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
                    EPSC_50 = two_pool(50, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC

                    r_squared_1hz = metrics.r2_score(EPSC_1hz[0:n_consider], EPSC_1[0:n_consider])
                    r_squared_10hz = metrics.r2_score(EPSC_10hz[0:n_consider], EPSC_10[0:n_consider])
                    r_squared_20hz = metrics.r2_score(EPSC_20hz[0:n_consider], EPSC_20[0:n_consider])
                    r_squared_50hz = metrics.r2_score(EPSC_50hz[0:n_consider], EPSC_50[0:n_consider])

                    avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4

                    if avg_r_squared > best[-1]:
                        best = [EPSC_1, EPSC_10, EPSC_20, EPSC_50, [p_fast, T_fast, T_slow, size_fast], [r_squared_1hz, r_squared_10hz, r_squared_20hz, r_squared_50hz], avg_r_squared]

    fig, axs = plt.subplots(4)

    l, = axs[0].plot(range(n_pulses), best[0])
    axs[0].scatter(range(n_pulses), EPSC_1hz)
    axs[0].set_xlim(0,n_pulses)
    axs[0].set_title('1 hz, R squared:' + str(best[-2][0]))
    axs[0].set_ylim(0,1)
    axs[0].set_ylabel('$EPSC/EPSC_{1}$')
    axs[0].set_xticks(range(n_pulses))

    m, = axs[1].plot(range(n_pulses), best[1])
    axs[1].scatter(range(n_pulses), EPSC_10hz)
    axs[1].set_xlim(0,n_pulses)
    axs[1].set_title('10 hz, R squared:' + str(best[-2][1]))
    axs[1].set_ylim(0,1)
    axs[1].set_ylabel('$EPSC/EPSC_{1}$')
    axs[1].set_xticks(range(n_pulses))

    n, = axs[2].plot(range(n_pulses), best[2])
    axs[2].scatter(range(n_pulses), EPSC_20hz)
    axs[2].set_xlim(0,n_pulses)
    axs[2].set_title('20 hz, R squared:' + str(best[-2][2]))
    axs[2].set_ylim(0,1)
    axs[2].set_ylabel('$EPSC/EPSC_{1}$')
    axs[2].set_xticks(range(n_pulses))

    o, = axs[3].plot(range(n_pulses), best[3])
    axs[3].scatter(range(n_pulses), EPSC_50hz)
    axs[3].set_xlim(0,n_pulses)
    axs[3].set_xlabel("Pulse #")
    axs[3].set_title('50 hz, R squared:' + str(best[-2][3]))
    axs[3].set_xticks(range(n_pulses))
    axs[3].set_ylim(0,1)
    axs[3].set_ylabel('$EPSC/EPSC_{1}$')
    fig.tight_layout()

    print('Best average R squared:', best[-1])
    print(best[-3])
    plt.show()
