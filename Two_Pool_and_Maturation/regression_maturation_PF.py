import numpy as np
import pandas as pd
import os
from TwoPoolAndMat import maturation
from matplotlib.widgets import Slider
from matplotlib import pyplot as plt
from sklearn import metrics

if __name__ == "__main__":

    data_1hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 1, nrows = 20)
    data_2hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 24, nrows = 20)
    data_10hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 47, nrows = 20)
    data_20hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 70, nrows = 20)
    data_50hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 93, nrows = 20)
    data_100hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 116, nrows = 20)

    EPSC_1hz = data_1hz.to_numpy()[:,0]
    stdev_1hz = data_1hz.to_numpy()[:,1]
    EPSC_2hz = data_2hz.to_numpy()[:,0]
    stdev_2hz = data_2hz.to_numpy()[:,1]
    EPSC_10hz = data_10hz.to_numpy()[:,0]
    stdev_10hz = data_10hz.to_numpy()[:,1]
    EPSC_20hz = data_20hz.to_numpy()[:,0]
    stdev_20hz = data_20hz.to_numpy()[:,1]
    EPSC_50hz = data_50hz.to_numpy()[:,0]
    stdev_50hz = data_50hz.to_numpy()[:,1]
    EPSC_100hz = data_100hz.to_numpy()[:,0]
    stdev_100hz = data_100hz.to_numpy()[:,1]

    n_sites = 1
    n_pulses = 20
    r = 20
    delta_t = 1000/r

    #release probabilities
    p_immature = .01
    p_mature = .05
    p_facilitated = .1

    #time constants
    T_refill = 160
    T_maturation = 2220
    T_facilitation = 160

    best = [-1e30]

    search_T_refill = np.linspace(100,300,10)
    search_T_maturation = np.linspace(1000,3000,10)
    search_p_mature = np.linspace(0.01,0.05,10)
    n_consider = 20

    for T_refill in search_T_refill:

        for T_maturation in search_T_maturation:

            for p_mature in search_p_mature:
                search_p_immature = np.linspace(0,p_mature,10)
                search_p_facilitated = np.linspace(p_mature, .5, 10)
                for p_immature in search_p_immature:

                    for p_facilitated in search_p_facilitated:

                        EPSC_1 = maturation(1, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
                        EPSC_2 = maturation(2, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
                        EPSC_10 = maturation(10, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
                        EPSC_20 = maturation(20, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
                        EPSC_50 = maturation(50, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
                        EPSC_100 = maturation(100, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC

                        r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
                        r_squared_2hz = metrics.r2_score(EPSC_2hz, EPSC_2)
                        r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
                        r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
                        r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)
                        r_squared_100hz = metrics.r2_score(EPSC_100hz, EPSC_100)

                        avg_r_squared = (r_squared_1hz + r_squared_2hz + r_squared_10hz + r_squared_20hz + r_squared_50hz + r_squared_100hz)/6

                        if avg_r_squared > best[-1]:
                            best = [EPSC_1, EPSC_2, EPSC_10, EPSC_20, EPSC_50, EPSC_100, [T_refill, T_maturation, p_immature, p_facilitated], [r_squared_1hz, r_squared_2hz, r_squared_10hz, r_squared_20hz, r_squared_50hz, r_squared_100hz], avg_r_squared]

    fig, axs = plt.subplots(6)

    l, = axs[0].plot(range(n_pulses), best[0])
    axs[0].scatter(range(n_pulses), EPSC_1hz)
    axs[0].set_xlim(0,n_pulses)
    axs[0].set_title('1 hz, R squared:' + str(best[-2][0]))
    axs[0].set_ylim(0,3)
    axs[0].set_ylabel('$EPSC/EPSC_{1}$')
    axs[0].set_xticks(range(n_pulses))

    m, = axs[1].plot(range(n_pulses), best[1])
    axs[1].scatter(range(n_pulses), EPSC_2hz)
    axs[1].set_xlim(0,n_pulses)
    axs[1].set_title('2 hz, R squared:' + str(best[-2][1]))
    axs[1].set_ylim(0,3)
    axs[1].set_ylabel('$EPSC/EPSC_{1}$')
    axs[1].set_xticks(range(n_pulses))

    n, = axs[2].plot(range(n_pulses), best[2])
    axs[2].scatter(range(n_pulses), EPSC_10hz)
    axs[2].set_xlim(0,n_pulses)
    axs[2].set_title('10 hz, R squared:' + str(best[-2][2]))
    axs[2].set_ylim(0,3)
    axs[2].set_ylabel('$EPSC/EPSC_{1}$')
    axs[2].set_xticks(range(n_pulses))

    o, = axs[3].plot(range(n_pulses), best[3])
    axs[3].scatter(range(n_pulses), EPSC_20hz)
    axs[3].set_xlim(0,n_pulses)
    axs[3].set_title('20 hz, R squared:' + str(best[-2][3]))
    axs[3].set_xticks(range(n_pulses))
    axs[3].set_ylim(0,3)
    axs[3].set_ylabel('$EPSC/EPSC_{1}$')

    p, = axs[4].plot(range(n_pulses), best[4])
    axs[4].scatter(range(n_pulses), EPSC_50hz)
    axs[4].set_xlim(0,n_pulses)
    axs[4].set_title('50 hz, R squared:' + str(best[-2][4]))
    axs[4].set_xticks(range(n_pulses))
    axs[4].set_ylim(0,3)
    axs[4].set_ylabel('$EPSC/EPSC_{1}$')

    q, = axs[5].plot(range(n_pulses), best[5])
    axs[5].scatter(range(n_pulses), EPSC_100hz)
    axs[5].set_xlim(0,n_pulses)
    axs[5].set_xlabel("Pulse #")
    axs[5].set_title('100 hz, R squared:' + str(best[-2][5]))
    axs[5].set_xticks(range(n_pulses))
    axs[5].set_ylim(0,3)
    axs[5].set_ylabel('$EPSC/EPSC_{1}$')
    fig.tight_layout()
    print('Best average R squared:', best[-1])
    print(best[-3])
    plt.show()
