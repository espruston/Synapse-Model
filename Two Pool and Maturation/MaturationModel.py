import numpy as np
import pandas as pd
import os
from TwoPoolAndMat import maturation
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

    best = [-1e30]
    search_T_maturation = np.linspace(1,1000,10)
    search_p_immature = np.linspace(0,0.9,10)
    search_p_mature = np.linspace(.6,.9,10)

    n_consider = 20
    r_squareds = []

    output = maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)
    print(metrics.r2_score(EPSC_20hz[0:n_consider], output.EPSC[0:n_consider]))

    # for T_maturation in search_T_maturation:
    #     for p_mature in search_p_mature:
    #         for p_immature in search_p_immature:
    #
    #             output = maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)
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

    l, = plt.plot(range(n_pulses), output.EPSC)
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

        output = maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)

        print(metrics.r2_score(EPSC_20hz[0:n_consider], output.EPSC[0:n_consider]))
        l.set_ydata(output.EPSC)
        plt.ylim(0,max(output.EPSC))
        fig.canvas.draw_idle()

    spim.on_changed(update)
    spma.on_changed(update)
    spfa.on_changed(update)
    sT_re.on_changed(update)
    sT_ma.on_changed(update)
    sT_fa.on_changed(update)

    plt.show()
