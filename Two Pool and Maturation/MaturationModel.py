import numpy as np
import pandas as pd
import os
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from math import e

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

    n_sites = 1
    n_pulses = 20
    r = 20
    delta_t = 1000/r

    #release probabilities
    p_immature = .2
    p_mature = .6
    p_facilitated = 1

    #time constants
    T_refill = 12
    T_maturation = 250
    T_facilitation = 2000

    # transition_mat = [1/T_refill, 0, 0, 0, -1/T_refill, 1/T_maturation+1/T_facilitation, 0, 0, 0, -1/T_maturation, 1/T_facilitation, 0, 0, -1/T_facilitation, -1/T_facilitation, 0] #markov transition matrix for 4 state system with facil. from I and M sites

    # transition_mat = [e**(-1*delta_t/T_refill), 0, 0, 0, 1-e**(-1*delta_t/T_refill), e**(-1*delta_t/T_maturation)*e**(-1*delta_t/T_facilitation), 0, 0, 0, 1-e**(-1*delta_t/T_maturation), e**(-1*delta_t/T_facilitation), 0, 0, 1-e**(-1*delta_t/T_facilitation), 1-e**(-1*delta_t/T_facilitation), 0]

    #transition_mat = np.reshape(transition_mat, (4,4))

    #release list preallocation
    n_release = []

    state = np.asarray([0, 0, n_sites, 0], dtype = 'float64') #initial state (empty, immature, mature, facil. sites)
    state = np.reshape(state, (4,1))
    EPSC = []

    best = [-1e30]
    search_T_maturation = np.linspace(1,1000,10)
    search_p_immature = np.linspace(0,0.9,10)
    search_p_mature = np.linspace(.6,.9,10)

    n_consider = 20
    r_squareds = []

    for i in range(n_pulses):

        n_release.append([state[1]*p_immature, state[2]*p_mature, state[3]*p_facilitated]) #immature, mature, and facilitated release
        EPSC.append(sum(n_release[-1]))

        state[0] += sum(n_release[-1]) #add released vesicles to empty sites
        state[1:] -= n_release[-1] #account for released vesicles in state change

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

    print(1 - sum((np.asarray(EPSC[0:n_consider])/EPSC[0] - np.average(EPSC_20hz[0:n_consider]))**2)/sum((EPSC_20hz[0:n_consider] - np.average(EPSC_20hz[0:n_consider]))**2))

    # for T_maturation in search_T_maturation:
    #     for p_mature in search_p_mature:
    #         for p_immature in search_p_immature:
    #
    #             n_release = []
    #
    #             state = np.asarray([0, 0, n_sites, 0], dtype = 'float64') #initial state (empty, immature, mature, facil. sites)
    #             state = np.reshape(state, (4,1))
    #
    #             EPSC = []
    #
    #             for i in range(n_pulses):
    #
    #                 n_release.append([state[1]*p_immature, state[2]*p_mature, state[3]*p_facilitated]) #immature, mature, and facilitated release
    #                 EPSC.append(sum(n_release[-1]))
    #
    #                 state[0] += sum(n_release[-1]) #add released vesicles to empty sites
    #                 state[1:] -= n_release[-1] #account for released vesicles in state change
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
    #                 r_squared = 1 - sum((np.asarray(EPSC[0:n_consider])/EPSC[0] - np.average(EPSC_20hz[0:n_consider]))**2)/sum((EPSC_20hz[0:n_consider] - np.average(EPSC_20hz[0:n_consider]))**2)
    #
    #             if r_squared > best[-1]:
    #                 best = [EPSC, T_maturation, T_refill, p_mature, p_immature, r_squared]
    #
    #                 r_squareds.append(r_squared)
    #
    # pathset = os.path.expanduser(r"~/Dropbox/Work/Jackman Lab/Modeling/200702_1_20hz_0FMM")
    # np.savez(pathset, best[0], r_squareds, best[-5:-1], best[-1])
    #
    # stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
    # max_time = int(1000/r*(n_pulses+3))
    # times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

    # alpha_1 = (times*e/2)*e**(-1*times/2) #reference alpha function
    # alpha = np.zeros(max_time) #functional alpha function
    #
    # for i in range(n_pulses):
    #     stimulus = stimulus_times[i]
    #     alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (best[0][i] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them
    # #
    # fig = plt.plot(times, -1*alpha/alpha[stimulus_times[0]])
    # plt.errorbar(stimulus_times, -1*EPSC_20hz, yerr = stdev_20hz, fmt = '.r', ecolor = 'black', elinewidth = 10 )
    #print(best)
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    l, = plt.plot(range(n_pulses), np.asarray(EPSC)/EPSC[0])
    plt.xlim(0,n_pulses+3)
    plt.ylim(0,1)
    plt.scatter(range(n_pulses), EPSC_20hz)

    axcolor = 'lightgoldenrodyellow'
    axpim = plt.axes([0.25, 0.025, 0.65, 0.01], facecolor=axcolor)
    axpma= plt.axes([0.25, 0.05, 0.65, 0.01], facecolor=axcolor)
    axpfa= plt.axes([0.25, 0.075, 0.65, 0.01], facecolor=axcolor)
    axT_re= plt.axes([0.25, 0.1, 0.65, 0.01], facecolor=axcolor)
    axT_ma= plt.axes([0.25, 0.125, 0.65, 0.01], facecolor=axcolor)
    axT_fa= plt.axes([0.25, 0.15, 0.65, 0.01], facecolor=axcolor)

    spim = Slider(axpim, 'p_immature', 0, 1, valinit=p_immature, valstep=.01)
    spma = Slider(axpma, 'p_mature', .1, 1, valinit=p_mature, valstep=.01)
    spfa = Slider(axpfa, 'p_facilitated', .1, 1, valinit=p_facilitated, valstep=.01)
    sT_re = Slider(axT_re, 'T_refill', 1, 2000, valinit=T_refill, valstep=20)
    sT_ma = Slider(axT_ma, 'T_maturation', 1, 2000, valinit=T_maturation, valstep=20)
    sT_fa = Slider(axT_fa, 'T_facilitation', 1, 2000, valinit=T_facilitation, valstep=20)

    def update(val):
        p_immature = spim.val
        p_mature = spma.val
        p_facilitated = spfa.val
        T_refill = sT_re.val
        T_maturation = sT_ma.val
        T_facilitation = sT_fa.val

        n_release = []

        state = np.asarray([0, 0, n_sites, 0], dtype = 'float64') #initial state (empty, immature, mature, facil. sites)
        state = np.reshape(state, (4,1))
        EPSC = []

        for i in range(n_pulses):
            n_release.append([state[1]*p_immature, state[2]*p_mature, state[3]*p_facilitated]) #immature, mature, and facilitated release
            EPSC.append(sum(n_release[-1]))

            state[0] += sum(n_release[-1]) #add released vesicles to empty sites
            state[1:] -= n_release[-1] #account for released vesicles in state change

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

        print(1 - sum((np.asarray(EPSC[0:n_consider])/EPSC[0] - np.average(EPSC_20hz[0:n_consider]))**2)/sum((EPSC_20hz[0:n_consider] - np.average(EPSC_20hz[0:n_consider]))**2))
        l.set_ydata(np.asarray(EPSC)/EPSC[0])
        plt.ylim(0,max(np.asarray(EPSC)/EPSC[0]))
        fig.canvas.draw_idle()

    spim.on_changed(update)
    spma.on_changed(update)
    spfa.on_changed(update)
    sT_re.on_changed(update)
    sT_ma.on_changed(update)
    sT_fa.on_changed(update)

    plt.show()
