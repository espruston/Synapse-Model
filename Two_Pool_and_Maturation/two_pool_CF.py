import numpy as np
import pandas as pd
from TwoPoolAndMat import two_pool
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

    r = 1

    m = .89
    size_fast = m
    size_slow = 1-size_fast

    p_fast = 0.56
    p_slow = 1-.198

    T_fast = .2
    T_slow = 7

    F = 1
    K_F = 1.5
    delta_F = 0
    T_F = 0.07
    BG_F = 1e-5
    sat_F = K_F
    facil = BG_F

    n_pulses = 20
    n_consider = 20

    output = two_pool(r, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    l, = plt.plot(range(n_pulses), output.EPSC, label = "EPSC")
    m, = plt.plot(range(n_pulses), output.fastEPSC, label = "Fast pool")
    n, = plt.plot(range(n_pulses), output.slowEPSC, label = "Slow pool")
    q, = plt.plot(range(n_pulses), EPSC_1hz, 'b.')
    plt.xlim(0,n_pulses)
    plt.ylim(0,1)

    axcolor = 'lightgoldenrodyellow'
    axpf = plt.axes([0.25, 0.025, 0.65, 0.01], facecolor=axcolor)
    axps= plt.axes([0.25, 0.05, 0.65, 0.01], facecolor=axcolor)
    axT_f= plt.axes([0.25, 0.075, 0.65, 0.01], facecolor=axcolor)
    axT_s= plt.axes([0.25, 0.1, 0.65, 0.01], facecolor=axcolor)
    axini_f= plt.axes([0.25, 0.125, 0.65, 0.01], facecolor=axcolor)
    axT_fa= plt.axes([0.25, 0.15, 0.65, 0.01], facecolor=axcolor)

    spf = Slider(axpf, 'p_fast', 0, 1, valinit=p_fast, valstep=.01)
    sps = Slider(axps, 'p_slow', .1, 1, valinit=p_slow, valstep=.01)
    sT_f = Slider(axT_f, 'T_fast', 0.01, 10, valinit=T_fast, valstep=0.01)
    sT_s = Slider(axT_s, 'T_slow', 0.01, 10, valinit=T_slow, valstep=0.01)
    sini_f = Slider(axini_f, 'size_fast', 0, 1, valinit=size_fast, valstep=.01)
    sT_fa = Slider(axT_fa, 'T_facilitation', 0.01, 10, valinit=T_F, valstep=0.01)

    rax = plt.axes([0.025, 0.5, 0.15, 0.20], facecolor=axcolor)
    radio = RadioButtons(rax, ('1hz', '10hz', '20hz', '50hz'), active=0)
    plt.title("Now showing data for:")

    def update(val):
        p_fast = spf.val
        p_slow = sps.val
        T_fast = sT_f.val
        T_slow = sT_s.val
        size_fast = sini_f.val
        size_slow = 1-size_fast
        T_F = sT_fa.val

        r = rs[radio.value_selected]

        output = two_pool(r, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F)

        EPSC_1 = two_pool(1, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
        EPSC_10 = two_pool(10, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
        EPSC_20 = two_pool(20, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC
        EPSC_50 = two_pool(50, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F).EPSC

        r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
        r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
        r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
        r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)

        avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4

        print(avg_r_squared)
        l.set_ydata(output.EPSC)
        m.set_ydata(output.fastEPSC)
        n.set_ydata(output.slowEPSC)
        q.set_ydata(EPSClabels[radio.value_selected])
        fig.canvas.draw_idle()

    spf.on_changed(update)
    sps.on_changed(update)
    sT_f.on_changed(update)
    sT_s.on_changed(update)
    sini_f.on_changed(update)
    sT_fa.on_changed(update)
    radio.on_clicked(update)

    plt.xlim(0,n_pulses)

    plt.show()
