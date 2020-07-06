import numpy as np
import pandas as pd
import os
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib import pyplot as plt
from sklearn import metrics
from math import e

if __name__ == "__main__":

    METHODS = {'regular': regular_train,
               'poisson': poisson_train,
               'RK1': DittmanRK1,
               'RK45': DittmanRK45}

    method = 'regular'
    method = METHODS[method]

    #arguments passed to simulator
    n_pulses = 20 #only needed for regular train
    r = 50 #frequency in Hz, only needed for regular train
    if method == METHODS['regular']:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int) #evenly spaced stimuli at r Hz starting at t = 1000/r ms
    else:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
        #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector [0,1,2,...,max_time]
    max_time = int(1000/r*(n_pulses+3)) #msec, used in plotting
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

    #climbing fiber analysis
    #file in /home/evan/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx

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

    n_consider = 20

    avg_1hz = np.average(EPSC_1hz[0:n_consider])
    avg_10hz = np.average(EPSC_10hz[0:n_consider])
    avg_20hz = np.average(EPSC_20hz[0:n_consider])
    avg_50hz = np.average(EPSC_50hz[0:n_consider])

    ss_1hz = sum((EPSC_1hz[0:n_consider] - avg_1hz)**2)
    ss_10hz = sum((EPSC_10hz[0:n_consider] - avg_10hz)**2)
    ss_20hz = sum((EPSC_20hz[0:n_consider] - avg_20hz)**2)
    ss_50hz = sum((EPSC_50hz[0:n_consider] - avg_50hz)**2)

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
    max_time4 = int(1000/r4*(n_pulses+3))

    def find_r_squared(args):

        EPSC_1 = method(stimulus_times1, max_time1, N_T, 0, F_1, args[0], args[1], args[2], args[3], args[4], 0, 1, T_E).EPSC
        EPSC_10 = method(stimulus_times2, max_time2, N_T, 0, F_1, args[0], args[1], args[2], args[3], args[4], 0, 1, T_E).EPSC
        EPSC_20 = method(stimulus_times3, max_time3, N_T, 0, F_1, args[0], args[1], args[2], args[3], args[4], 0, 1, T_E).EPSC
        EPSC_50 = method(stimulus_times4, max_time4, N_T, 0, F_1, args[0], args[1], args[2], args[3], args[4], 0, 1, T_E).EPSC

        r_squared_1hz = metrics.r2_score(EPSC_1hz, EPSC_1)
        r_squared_10hz = metrics.r2_score(EPSC_10hz, EPSC_10)
        r_squared_20hz = metrics.r2_score(EPSC_20hz, EPSC_20)
        r_squared_50hz = metrics.r2_score(EPSC_50hz, EPSC_50)

        avg_r_squared = (r_squared_1hz + r_squared_10hz + r_squared_20hz + r_squared_50hz)/4

        vec = [EPSC_1, EPSC_10, EPSC_20, EPSC_50, args, [r_squared_1hz, r_squared_10hz, r_squared_20hz, r_squared_50hz], avg_r_squared]

        return(vec)

    """brute force search"""
    best = [-1e309]
    search_T_D = np.linspace(1,350,100)
    search_k_0 = np.linspace(.1, 10.1, 100)
    search_k_max = np.linspace(0.1, 20, 100)

    for T_D in search_T_D:
        for k_0 in search_k_0:
            for k_max in search_k_max:

                output = find_r_squared([T_F, T_D, K_D, k_0, k_max])
                if output[-1] > best[-1]:
                    best = output.copy()

    """linear descent search"""
    # args = [T_F, T_D, K_D, k_0, k_max]
    # args_plus = [T_F, T_D, K_D, k_0, k_max]
    # args_minus = [T_F, T_D, K_D, k_0, k_max]
    # args_min = [1, 1, .1, .1, .1]
    # args_max = [1000, 1000, 5, 20, 100]
    # deltas = [1, 1, .1, .1, .1]
    # best = find_r_squared(args)
    #
    # c = 0
    # while c < 10000:
    #     c_before = c
    #     i = 0
    #     for i in range(len(args)):
    #         while args_plus[i] + deltas[i] < args_max[i] and args_minus[i] - deltas[i] > args_min[i]:
    #
    #             args_plus[i] += deltas[i]
    #             args_minus[i] -= deltas[i]
    #
    #             plus = find_r_squared(args_plus)
    #             minus = find_r_squared(args_minus)
    #
    #             c += 1
    #
    #             if plus[-1] > best[-1]:
    #                 if minus[-1] > plus[-1]: #both are improvements, but minus is bigger than plus
    #                     best = minus.copy()
    #                     args = args_minus.copy()
    #                     args_plus = args.copy()
    #                 else: #plus is a bigger improvement than minus
    #                     best = plus.copy()
    #                     args = args_plus.copy()
    #                     args_minus = args.copy()
    #             elif minus[-1] > best[-1]: #minus is an improvement but plus isnt
    #                 best = minus.copy()
    #                 args = args_minus.copy()
    #                 args_plus = args.copy()
    #             else: #neither plus or minus improves, so move to next variable
    #                 break
    #             if c%10 == 0:
    #                 print(best[-1])
    #
    #     if c == c_before:
    #         print("No improvements available")
    #         break

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

    print('R squared:', best[-1])
    print(best[-3])
    plt.show()
