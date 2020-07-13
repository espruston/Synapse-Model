import numpy as np
import pandas as pd
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib import pyplot as plt
from sklearn import metrics
from math import e

if __name__ == "__main__":

    data_WT = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Raw data climbing fiber bursts DW.xlsx', usecols = 'A,B,C', skiprows = 93, nrows = 11)

    ITIs = data_WT.to_numpy()[:,0]
    EPSCs2_WT = data_WT.to_numpy()[:,1]
    stdev_WT = data_WT.to_numpy()[:,2]

    METHODS = {'regular': regular_train,
               'poisson': poisson_train,
               'RK1': DittmanRK1,
               'RK45': DittmanRK45}

    method = 'RK1'
    method = METHODS[method]

    r = 50

    N_T = 1 #number of total release sites
    roh = 0 #EPSC2/EPSC1
    F_1 = 1-.198 #initial release probability
    T_F = 100 #decay constant for CaX_F
    T_D = 50 #decay constant for CaX_D
    K_D = 2 #affinity of CaX_D for site
    k_0 = 1 #initial recovery rate
    k_max = 20 #maximum recovery rate
    #K_F = 1 #affinity of CaX_F for site
    delta_F = 0 #amount by which CaX_F increases as a result of stimulus
    delta_D = 1 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    n_pulses = 20

    EPSCs2_reg = []

    for t in ITIs:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
        stimulus_times.append(simutulus_times[-1]+t)
        max_time = int(1000/r*(n_pulses+3))+t

        output = method(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

        EPSCs2_reg.append(output.EPSC[-1])

    plt.errorbar(ITIs, EPSCs2_WT, yerr = stdev_WT, label = "data")
    plt.plot(ITIs, EPSCs2_reg, label = "regression")
    print(metrics.r2_score(EPSCs2_WT, EPSCs2_reg))

    # print('Best R squared:', best[-1])
    # print(best[-2])
    plt.legend()
    plt.show()
