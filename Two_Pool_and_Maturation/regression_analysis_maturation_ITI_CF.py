import numpy as np
import pandas as pd
from TwoPoolAndMat import maturation
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


    n_sites = 1
    n_pulses = 20
    r = 50

    #release probabilities
    p_immature = .624
    p_mature = 1-.198
    p_facilitated = 1

    #time constants
    T_refill = 160
    T_maturation = 2200
    T_facilitation = 1e30

    n_consider = 20

    best = [-1e30]
    search_T_maturation = np.linspace(1000,3000,10)
    search_T_refill = np.linspace(100,300,10)
    search_p_immature = np.linspace(0,p_mature,10)

    n_consider = 20
    EPSCs2_reg = []
    first_burst = maturation(50, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)

    for t in ITIs:
        state = first_burst.final_state
        delta_t = t - 1000/r #1 ISI has already been simulated for the final state

        n_E_I = state[0]*(1-e**(-1*delta_t/T_refill)) #start in E end in I
        n_E_I_F = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_facilitation)) #start in E end in F after passing through I
        n_E_I_M = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_maturation))
        n_E_I_M_F = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_maturation))*(1-e**(-1*delta_t/T_facilitation))
        n_I_M = state[1]*(1-e**(-1*delta_t/T_maturation))
        n_I_F = state[1]*(1-e**(-1*delta_t/T_facilitation))
        n_I_M_F =  state[1]*(1-e**(-1*delta_t/T_maturation))*(1-e**(-1*delta_t/T_facilitation))
        n_M_F = state[2]*(1-e**(-1*delta_t/T_facilitation))

        state[0] = state[0]*e**(-1*delta_t/T_refill) #after vesicles are released, additional empty sites refill with characteristic time T_refill over delta_t, reducing the number of empty sites
        state[1] += n_E_I - n_I_M - n_E_I_M - n_I_F - n_E_I_F  #the pool begins to refill with immature vesicles, over delta_t, some facilitate, some mature, and a small number do both, some of the old vesicles mature and some facilitate

        state[2] += n_I_M + n_E_I_M - n_M_F - n_I_M_F - n_E_I_M_F #some of the mature pool is replenished by the immature pool that remains after the pulse and some is replenished by the new vesicles refilling pool first, some mature vesicles are lost to facilitation

        state[3] += n_E_I_F + n_E_I_M_F + n_I_F + n_I_M_F + n_M_F #vesicles become facilitated by 5 paths, 2 starting empty, 2 starting immature, and one starting mature

        EPSCs2_reg.append((state[1]*p_immature + state[2]*p_mature + state[3]*p_facilitated)/first_burst.norm_val)
    # for T_refill in search_T_refill:
    #     for T_maturation in search_T_maturation:
    #         for p_immature in search_p_immature:
    #             EPSCs2_reg = []
    #             first_burst = maturation(50, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)
    #
    #             for t in ITIs:
    #                 state = first_burst.final_state
    #                 delta_t = t - 1000/r #1 ISI has already been simulated for the final state
    #
    #                 n_E_I = state[0]*(1-e**(-1*delta_t/T_refill)) #start in E end in I
    #                 n_E_I_F = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_facilitation)) #start in E end in F after passing through I
    #                 n_E_I_M = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_maturation))
    #                 n_E_I_M_F = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_maturation))*(1-e**(-1*delta_t/T_facilitation))
    #                 n_I_M = state[1]*(1-e**(-1*delta_t/T_maturation))
    #                 n_I_F = state[1]*(1-e**(-1*delta_t/T_facilitation))
    #                 n_I_M_F =  state[1]*(1-e**(-1*delta_t/T_maturation))*(1-e**(-1*delta_t/T_facilitation))
    #                 n_M_F = state[2]*(1-e**(-1*delta_t/T_facilitation))
    #
    #                 state[0] = state[0]*e**(-1*delta_t/T_refill) #after vesicles are released, additional empty sites refill with characteristic time T_refill over delta_t, reducing the number of empty sites
    #                 state[1] += n_E_I - n_I_M - n_E_I_M - n_I_F - n_E_I_F  #the pool begins to refill with immature vesicles, over delta_t, some facilitate, some mature, and a small number do both, some of the old vesicles mature and some facilitate
    #
    #                 state[2] += n_I_M + n_E_I_M - n_M_F - n_I_M_F - n_E_I_M_F #some of the mature pool is replenished by the immature pool that remains after the pulse and some is replenished by the new vesicles refilling pool first, some mature vesicles are lost to facilitation
    #
    #                 state[3] += n_E_I_F + n_E_I_M_F + n_I_F + n_I_M_F + n_M_F #vesicles become facilitated by 5 paths, 2 starting empty, 2 starting immature, and one starting mature
    #
    #                 EPSCs2_reg.append((state[1]*p_immature + state[2]*p_mature + state[3]*p_facilitated)/first_burst.norm_val)
    #
    #             r_squared = metrics.r2_score(EPSCs2_WT, EPSCs2_reg)
    #
    #             if r_squared > best[-1]:
    #                 best = [EPSCs2_reg, [T_refill, T_maturation, p_immature], r_squared]
    #
    plt.errorbar(ITIs, EPSCs2_WT, yerr = stdev_WT, label = "data")
    plt.plot(ITIs, EPSCs2_reg, label = "regression")
    print(metrics.r2_score(EPSCs2_WT, EPSCs2_reg))

    # print('Best R squared:', best[-1])
    # print(best[-2])
    plt.legend()
    plt.show()
