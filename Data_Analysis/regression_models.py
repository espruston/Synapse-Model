import pandas as pd
import numpy as np
from math import log, e
from matplotlib import pyplot as plt
from sklearn import metrics


if __name__ == "__main__":

    #import raw data to pandas dataframes
    CF_data_1hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 1, nrows = 20)
    CF_data_10hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 23, nrows = 20)
    CF_data_20hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 45, nrows = 20)
    CF_data_50hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 67, nrows = 20)

    PF_data_1hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 1, nrows = 20)
    PF_data_10hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 47, nrows = 20)
    PF_data_20hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 70, nrows = 20)
    PF_data_50hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 93, nrows = 20)

    #define data series as numpy vectors
    CF_WT_1hz_EPSCs = CF_data_1hz.to_numpy()[:,0]
    CF_WT_1hz_SEMs = CF_data_1hz.to_numpy()[:,1]
    CF_KO_1hz_EPSCs = CF_data_1hz.to_numpy()[:,2]
    CF_KO_1hz_SEMs = CF_data_1hz.to_numpy()[:,3]
    CF_WTCaBlock_1hz_EPSCs = CF_data_1hz.to_numpy()[:,4]
    CF_WTCaBlock_1hz_SEMs = CF_data_1hz.to_numpy()[:,5]
    CF_KOCaBlock_1hz_EPSCs = CF_data_1hz.to_numpy()[:,6]
    CF_KOCaBlock_1hz_SEMs = CF_data_1hz.to_numpy()[:,7]

    CF_WT_10hz_EPSCs = CF_data_10hz.to_numpy()[:,0]
    CF_WT_10hz_SEMs = CF_data_10hz.to_numpy()[:,1]
    CF_KO_10hz_EPSCs = CF_data_10hz.to_numpy()[:,2]
    CF_KO_10hz_SEMs = CF_data_10hz.to_numpy()[:,3]
    CF_WTCaBlock_10hz_EPSCs = CF_data_10hz.to_numpy()[:,4]
    CF_WTCaBlock_10hz_SEMs = CF_data_10hz.to_numpy()[:,5]
    CF_KOCaBlock_10hz_EPSCs = CF_data_10hz.to_numpy()[:,6]
    CF_KOCaBlock_10hz_SEMs = CF_data_10hz.to_numpy()[:,7]

    CF_WT_20hz_EPSCs = CF_data_20hz.to_numpy()[:,0]
    CF_WT_20hz_SEMs = CF_data_20hz.to_numpy()[:,1]
    CF_KO_20hz_EPSCs = CF_data_20hz.to_numpy()[:,2]
    CF_KO_20hz_SEMs = CF_data_20hz.to_numpy()[:,3]
    CF_WTCaBlock_20hz_EPSCs = CF_data_20hz.to_numpy()[:,4]
    CF_WTCaBlock_20hz_SEMs = CF_data_20hz.to_numpy()[:,5]
    CF_KOCaBlock_20hz_EPSCs = CF_data_20hz.to_numpy()[:,6]
    CF_KOCaBlock_20hz_SEMs = CF_data_20hz.to_numpy()[:,7]

    CF_WT_50hz_EPSCs = CF_data_50hz.to_numpy()[:,0]
    CF_WT_50hz_SEMs = CF_data_50hz.to_numpy()[:,1]
    CF_KO_50hz_EPSCs = CF_data_50hz.to_numpy()[:,2]
    CF_KO_50hz_SEMs = CF_data_50hz.to_numpy()[:,3]
    CF_WTCaBlock_50hz_EPSCs = CF_data_50hz.to_numpy()[:,4]
    CF_WTCaBlock_50hz_SEMs = CF_data_50hz.to_numpy()[:,5]
    CF_KOCaBlock_50hz_EPSCs = CF_data_50hz.to_numpy()[:,6]
    CF_KOCaBlock_50hz_SEMs = CF_data_50hz.to_numpy()[:,7]

    PF_WT_1hz_EPSCs = PF_data_1hz.to_numpy()[:,0]
    PF_WT_1hz_SEMs = PF_data_1hz.to_numpy()[:,1]
    PF_KO_1hz_EPSCs = PF_data_1hz.to_numpy()[:,2]
    PF_KO_1hz_SEMs = PF_data_1hz.to_numpy()[:,3]
    PF_WTCaBlock_1hz_EPSCs = PF_data_1hz.to_numpy()[:,4]
    PF_WTCaBlock_1hz_SEMs = PF_data_1hz.to_numpy()[:,5]
    PF_KOCaBlock_1hz_EPSCs = PF_data_1hz.to_numpy()[:,6]
    PF_KOCaBlock_1hz_SEMs = PF_data_1hz.to_numpy()[:,7]

    PF_WT_10hz_EPSCs = PF_data_10hz.to_numpy()[:,0]
    PF_WT_10hz_SEMs = PF_data_10hz.to_numpy()[:,1]
    PF_KO_10hz_EPSCs = PF_data_10hz.to_numpy()[:,2]
    PF_KO_10hz_SEMs = PF_data_10hz.to_numpy()[:,3]
    PF_WTCaBlock_10hz_EPSCs = PF_data_10hz.to_numpy()[:,4]
    PF_WTCaBlock_10hz_SEMs = PF_data_10hz.to_numpy()[:,5]
    PF_KOCaBlock_10hz_EPSCs = PF_data_10hz.to_numpy()[:,6]
    PF_KOCaBlock_10hz_SEMs = PF_data_10hz.to_numpy()[:,7]

    PF_WT_20hz_EPSCs = PF_data_20hz.to_numpy()[:,0]
    PF_WT_20hz_SEMs = PF_data_20hz.to_numpy()[:,1]
    PF_KO_20hz_EPSCs = PF_data_20hz.to_numpy()[:,2]
    PF_KO_20hz_SEMs = PF_data_20hz.to_numpy()[:,3]
    PF_WTCaBlock_20hz_EPSCs = PF_data_20hz.to_numpy()[:,4]
    PF_WTCaBlock_20hz_SEMs = PF_data_20hz.to_numpy()[:,5]
    PF_KOCaBlock_20hz_EPSCs = PF_data_20hz.to_numpy()[:,6]
    PF_KOCaBlock_20hz_SEMs = PF_data_20hz.to_numpy()[:,7]

    PF_WT_50hz_EPSCs = PF_data_50hz.to_numpy()[:,0]
    PF_WT_50hz_SEMs = PF_data_50hz.to_numpy()[:,1]
    PF_KO_50hz_EPSCs = PF_data_50hz.to_numpy()[:,2]
    PF_KO_50hz_SEMs = PF_data_50hz.to_numpy()[:,3]
    PF_WTCaBlock_50hz_EPSCs = PF_data_50hz.to_numpy()[:,4]
    PF_WTCaBlock_50hz_SEMs = PF_data_50hz.to_numpy()[:,5]
    PF_KOCaBlock_50hz_EPSCs = PF_data_50hz.to_numpy()[:,6]
    PF_KOCaBlock_50hz_SEMs = PF_data_50hz.to_numpy()[:,7]

    ### CHANGE THESE TWO VALUES ###
    data = CF_WT_20hz_EPSCs
    hz = 20
    n_pulses = 20
    ### ----------------------- ###

    """simple exponential decay"""
    # x = np.linspace(0,19*1000/hz,20)
    #
    # lny = []
    # for i in range(len(data)):
    #     lny.append(log(data[i]))
    #
    # b = sum(x*lny)/sum(x**2)
    #
    # print('R squared:', metrics.r2_score(data, e**(b*x)))
    # print('Equation: y = e^(',b,'* x)')

    EPSC_1hz = CF_WT_1hz_EPSCs
    EPSC_10hz = CF_WT_10hz_EPSCs
    EPSC_20hz = CF_WT_20hz_EPSCs
    EPSC_50hz = CF_WT_50hz_EPSCs

    p_release = 1-.198
    T_refill = 250
    """constant release probability with constant refill"""
    def simple_refill(n_pulses, r, p_release, T_refill):
        delta_t = 1000/r
        RRP = 1
        EPSC = []
        for i in range(n_pulses):
            EPSC.append(p_release*RRP)
            RRP -= p_release*RRP

            RRP += (1-RRP)*(1-e**(-1*delta_t/T_refill))
        return(np.asarray(EPSC)/EPSC[0])

    """brute force search the consant refill case"""
    best = [-1e309]
    search_T_refill = np.linspace(1,1000,1000)
    search_p_release = np.linspace(.01,1,50)
    for T_refill in search_T_refill:
        for p_release in search_p_release:
            EPSC_1 = simple_refill(n_pulses, 1, p_release, T_refill)
            EPSC_10 = simple_refill(n_pulses, 10, p_release, T_refill)
            EPSC_20 = simple_refill(n_pulses, 20, p_release, T_refill)
            EPSC_50 = simple_refill(n_pulses, 50, p_release, T_refill)

            r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
            r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
            r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
            r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)

            avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4

            if avg_r_squared > best[-1]:
                best = [EPSC_1, EPSC_10, EPSC_20, EPSC_50, T_refill, p_release, [r_squared_1hz, r_squared_10hz, r_squared_20hz, r_squared_50hz], avg_r_squared]

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
    print(best[-4])

    # fig, ax1 = plt.subplots()
    # ax1.scatter((x+1000/hz)*hz/1000, data)
    # ax1.set_xlabel('Pulse #')
    # ax1.set_xlim(1,22)
    # ax1.set_xticks(np.linspace(1,20,20))
    # ax1.set_ylabel('$EPSC/EPSC_{1}$')
    # ax1.set_ylim(0,max(data)+.5)
    #
    # xlong = np.linspace(0,19*1000/hz+2*1000/hz,22)
    # ax2 = ax1.twiny()
    # ax2.plot(xlong,e**(b*xlong))
    # ax2.set_xlabel('t (ms)')
    # ax2.set_xlim(0,max(xlong))
    # ax2.set_xticks(xlong)
    #
    # plt.grid(axis = 'x')

    plt.show()
