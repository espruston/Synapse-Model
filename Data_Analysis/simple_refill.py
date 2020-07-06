import pandas as pd
import numpy as np
from math import log, e
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

    n_pulses = 20
    n_consider = 20
    #EPSC_1 = simple_refill(n_pulses, 1, p_release, T_refill)
    #EPSC_10 = simple_refill(n_pulses, 10, p_release, T_refill)
    EPSC_20 = simple_refill(n_pulses, 20, p_release, T_refill)
    #EPSC_50 = simple_refill(n_pulses, 50, p_release, T_refill)

    #r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
    #r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
    r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
    #r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)

    #avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4
    #print(avg_r_squared)
    print(r_squared_20hz)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    l, = plt.plot(range(n_pulses), EPSC_20)
    plt.xlim(0,n_pulses)
    plt.xticks(range(n_pulses))
    #plt.ylim(0,1)
    plt.scatter(range(n_pulses), EPSC_20hz)

    axcolor = 'lightgoldenrodyellow'
    axp = plt.axes([0.25, 0.025, 0.65, 0.01], facecolor=axcolor)
    axT = plt.axes([0.25, 0.05, 0.65, 0.01], facecolor=axcolor)

    sp = Slider(axp, 'p_release', 0, 1, valinit=p_release, valstep=.01)
    sT = Slider(axT, 'T_refill', 1, 1000, valinit=T_refill, valstep=20)

    def update(val):
        p_release = sp.val
        T_refill = sT.val
        r = 20

        EPSC = simple_refill(n_pulses, r, p_release, T_refill)
        print(metrics.r2_score(EPSC_20hz[0:n_consider], EPSC[0:n_consider]))
        l.set_ydata(EPSC)
        plt.ylim(0,max(EPSC))
        fig.canvas.draw()

    sp.on_changed(update)
    sT.on_changed(update)

    plt.show()
