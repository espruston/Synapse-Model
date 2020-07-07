import numpy as np
import pandas as pd
import os
from TwoPoolAndMat import maturation
from matplotlib.widgets import Slider, RadioButtons
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

    rs = {'1hz': 1, '10hz': 10,  '20hz': 20, '50hz': 50}

    EPSClabels = {'1hz': EPSC_1hz, '10hz': EPSC_10hz, '20hz': EPSC_20hz, '50hz': EPSC_50hz}

    n_sites = 1
    n_pulses = 20
    r = 1
    delta_t = 1000/r

    #release probabilities
    p_immature = .624
    p_mature = 1-.198
    p_facilitated = 1

    #time constants
    T_refill = 160
    T_maturation = 2200
    T_facilitation = 1e30

    n_consider = 20
    r_squareds = []

    output = maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    l, = plt.plot(range(n_pulses), output.EPSC, label = "EPSC")
    m, = plt.plot(range(n_pulses), output.EPSC_immature, label = "Immature")
    n, = plt.plot(range(n_pulses), output.EPSC_mature, label = "Mature")
    o, = plt.plot(range(n_pulses), output.EPSC_facilitated, label = "Facilitated")
    q, = plt.plot(range(n_pulses), EPSC_1hz, 'b.')
    plt.xlim(0,n_pulses)
    plt.ylim(0,1)


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
    sT_re = Slider(axT_re, 'T_refill', 1, 500, valinit=T_refill, valstep=20)
    sT_ma = Slider(axT_ma, 'T_maturation', 1000, 4000, valinit=T_maturation, valstep=20)
    sT_fa = Slider(axT_fa, 'T_facilitation', 1, 1e30, valinit=T_facilitation, valstep=20)

    rax = plt.axes([0.025, 0.5, 0.15, 0.20], facecolor=axcolor)
    radio = RadioButtons(rax, ('1hz', '10hz', '20hz', '50hz'), active=0)
    plt.title("Now showing data for:")

    def update(val):
        p_immature = spim.val
        p_mature = spma.val
        p_facilitated = spfa.val
        T_refill = sT_re.val
        T_maturation = sT_ma.val
        T_facilitation = sT_fa.val

        r = rs[radio.value_selected]

        EPSC_1 = maturation(1, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
        EPSC_10 = maturation(10, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
        EPSC_20 = maturation(20, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC
        EPSC_50 = maturation(50, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation).EPSC

        r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
        r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
        r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
        r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)

        avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4
        output = maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)

        print(avg_r_squared)
        l.set_ydata(output.EPSC)
        m.set_ydata(output.EPSC_immature)
        n.set_ydata(output.EPSC_mature)
        o.set_ydata(output.EPSC_facilitated)
        q.set_ydata(EPSClabels[radio.value_selected])
        fig.canvas.draw_idle()

    spim.on_changed(update)
    spma.on_changed(update)
    spfa.on_changed(update)
    sT_re.on_changed(update)
    sT_ma.on_changed(update)
    sT_fa.on_changed(update)
    radio.on_clicked(update)

    plt.show()
