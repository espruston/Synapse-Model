import numpy as np
import pandas as pd
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib.widgets import Slider
from matplotlib import pyplot as plt
from sklearn import metrics

if __name__ == "__main__":

    METHODS = {'regular': regular_train,
               'poisson': poisson_train,
               'RK1': DittmanRK1,
               'RK45': DittmanRK45}

    method = 'RK1'
    method = METHODS[method]

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

    #arguments passed to simulator
    n_pulses = 20 #only needed for regular train
    n_consider = 20
    r = 20 #frequency in Hz, only needed for regular train
    if method == METHODS['regular']:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int) #evenly spaced stimuli at r Hz starting at t = 1000/r ms
    else:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
        #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector [0,1,2,...,max_time]
    max_time = int(1000/r*(n_pulses+3)) #msec, used in plotting
    N_T = 1 #number of total release sites
    rho = 0 #EPSC2/EPSC1
    F_1 = 1-.198 #initial release probability
    T_F = 50 #decay constant for CaX_F
    T_D = 50 #decay constant for CaX_D
    K_D = 2 #affinity of CaX_D for site
    k_0 = 1.216 #initial recovery rate
    k_max = 20 #maximum recovery rate
    #K_F = 1 #affinity of CaX_F for site
    delta_F = 1 #amount by which CaX_F increases as a result of stimulus
    delta_D = 1 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    output = method(stimulus_times, max_time, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    print(metrics.r2_score(EPSC_20hz[0:n_consider], output.EPSC[0:n_consider]))

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    l, = plt.plot(range(n_pulses), output.EPSC)
    plt.xlim(0,n_pulses)
    plt.xticks(range(n_pulses))
    #plt.ylim(0,1)
    plt.scatter(range(n_pulses), EPSC_20hz)

    axcolor = 'lightgoldenrodyellow'
    axrho = plt.axes([0.25, 0.025, 0.65, 0.01], facecolor=axcolor)
    axF_1 = plt.axes([0.25, 0.05, 0.65, 0.01], facecolor=axcolor)
    axT_F = plt.axes([0.25, 0.075, 0.65, 0.01], facecolor=axcolor)
    axT_D = plt.axes([0.25, 0.1, 0.65, 0.01], facecolor=axcolor)
    axk_0 = plt.axes([0.25, 0.125, 0.65, 0.01], facecolor=axcolor)
    axk_max = plt.axes([0.25, 0.15, 0.65, 0.01], facecolor=axcolor)

    srho = Slider(axrho, 'rho', 0, 10, valinit=rho, valstep=.1)
    sF_1 = Slider(axF_1, 'F_1', 0, 1, valinit=F_1, valstep=.01)
    sT_F = Slider(axT_F, 'T_F', 1, 2000, valinit=T_F, valstep=20)
    sT_D = Slider(axT_D, 'T_D', 1, 2000, valinit=T_D, valstep=20)
    sk_0 = Slider(axk_0, 'k_0', 0, 20, valinit=k_0, valstep=.2)
    sk_max = Slider(axk_max, 'k_max', k_0, 100, valinit=k_max, valstep=.2)

    def update(val):
        rho = srho.val
        F_1 = sF_1.val
        T_F = sT_F.val
        T_D = sT_D.val
        k_0 = sk_0.val
        k_max = sk_max.val

        output = method(stimulus_times, max_time, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

        print(metrics.r2_score(EPSC_20hz[0:n_consider], output.EPSC[0:n_consider]))
        l.set_ydata(output.EPSC)
        plt.ylim(0,max(output.EPSC))
        fig.canvas.draw()

    srho.on_changed(update)
    sF_1.on_changed(update)
    sT_F.on_changed(update)
    sT_D.on_changed(update)
    sk_0.on_changed(update)
    sk_max.on_changed(update)
    plt.show()
