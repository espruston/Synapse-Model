import numpy as np
import pandas as pd
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib.widgets import Slider, RadioButtons
from matplotlib import pyplot as plt
from sklearn import metrics

if __name__ == "__main__":

    METHODS = {'regular': regular_train,
               'poisson': poisson_train,
               'RK1': DittmanRK1,
               'RK45': DittmanRK45}

    method = 'regular'
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
    r1 = 1
    stimulus_times1 = np.linspace(1000/r1,(1000/r1)*n_pulses,n_pulses, dtype = int)
    max_time1 = int(1000/r1*(n_pulses+3))
    r2 = 10
    stimulus_times2 = np.linspace(1000/r2,(1000/r2)*n_pulses,n_pulses, dtype = int)
    max_time2 = int(1000/r2*(n_pulses+3))
    r3 = 20
    stimulus_times3 = np.linspace(1000/r3,(1000/r3)*n_pulses,n_pulses, dtype = int)
    max_time3 = int(1000/r3*(n_pulses+3))
    r4 = 50
    stimulus_times4 = np.linspace(1000/r4,(1000/r4)*n_pulses,n_pulses, dtype = int)
    max_time4 = int(1000/r4*(n_pulses+3))    #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector [0,1,2,...,max_time]

    stimdict = {'1hz': stimulus_times1, '10hz': stimulus_times2,  '20hz': stimulus_times3, '50hz': stimulus_times4}
    max_timedict = {'1hz': max_time1, '10hz': max_time2,  '20hz': max_time3, '50hz': max_time4}

    EPSClabels = {'1hz': EPSC_1hz, '10hz': EPSC_10hz, '20hz': EPSC_20hz, '50hz': EPSC_50hz}

    N_T = 1 #number of total release sites
    rho = 0 #EPSC2/EPSC1
    F_1 = 1-.198 #initial release probability
    T_F = 50 #decay constant for CaX_F
    T_D = 96 #decay constant for CaX_D
    K_D = 2 #affinity of CaX_D for site
    k_0 = 0.302 #initial recovery rate
    k_max = 1.09 #maximum recovery rate
    #K_F = 1 #affinity of CaX_F for site
    delta_F = 0 #amount by which CaX_F increases as a result of stimulus
    delta_D = 1 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    output = method(stimulus_times1, max_time1, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    l, = plt.plot(range(n_pulses), output.EPSC)
    q, = plt.plot(range(n_pulses), EPSC_1hz, 'b.')
    plt.ylim(0,1)
    plt.xlim(0,n_pulses)
    plt.xticks(range(n_pulses))

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

    rax = plt.axes([0.025, 0.5, 0.15, 0.20], facecolor=axcolor)
    radio = RadioButtons(rax, ('1hz', '10hz', '20hz', '50hz'), active=0)
    plt.title("Now showing data for:")

    def update(val):
        rho = srho.val
        F_1 = sF_1.val
        T_F = sT_F.val
        T_D = sT_D.val
        k_0 = sk_0.val
        k_max = sk_max.val

        output = method(stimdict[radio.value_selected], max_timedict[radio.value_selected], N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

        EPSC_1 = method(stimulus_times1, max_time1, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E).EPSC
        EPSC_10 = method(stimulus_times2, max_time2, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E).EPSC
        EPSC_20 = method(stimulus_times3, max_time3, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E).EPSC
        EPSC_50 = method(stimulus_times4, max_time4, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E).EPSC

        r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
        r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
        r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
        r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)

        avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4
        print(avg_r_squared)

        l.set_ydata(output.EPSC)
        q.set_ydata(EPSClabels[radio.value_selected])
        fig.canvas.draw()

    srho.on_changed(update)
    sF_1.on_changed(update)
    sT_F.on_changed(update)
    sT_D.on_changed(update)
    sk_0.on_changed(update)
    sk_max.on_changed(update)
    radio.on_clicked(update)
    plt.show()
