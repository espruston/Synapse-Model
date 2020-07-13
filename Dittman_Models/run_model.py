import numpy as np
import pandas as pd
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
    r = 20 #frequency in Hz, only needed for regular train
    if method == METHODS['regular']:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int) #evenly spaced stimuli at r Hz starting at t = 1000/r ms
    else:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
        #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector [0,1,2,...,max_time]
    max_time = int(1000/r*(n_pulses+3)) #msec, used in plotting
    N_T = 1 #number of total release sites
    roh = 0 #EPSC2/EPSC1
    F_1 = 1-.198 #initial release probability
    T_F = 50 #decay constant for CaX_F
    T_D = 96.2 #decay constant for CaX_D
    K_D = 2 #affinity of CaX_D for site
    k_0 = 0.302 #initial recovery rate
    k_max = 1.09 #maximum recovery rate
    #K_F = 1 #affinity of CaX_F for site
    delta_F = 0 #amount by which CaX_F increases as a result of stimulus
    delta_D = 1 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    """TEST ONE METHOD ON ONE SET OF INPUTS"""
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

    #output = method(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    # print(output.EPSC)
    # print(EPSC_20hz)
    # print(metrics.r2_score(EPSC_20hz, output.EPSC))
    # fig = plt.plot(range(n_pulses), output.EPSC)
    # plt.plot(range(n_pulses), EPSC_20hz, 'b.')
    # plt.show()
    """TEST TWO METHODS"""
    output1 = regular_train(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    output2 = DittmanRK1(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    # print(np.asarray(output1.F)*np.asarray(output1.D)/.802)
    print(output1.EPSC)
    # print(output2.F[stimulus_times-1]*output2.D[stimulus_times-1]/.802)
    print(output2.EPSC)
    #print(output2.F[stimulus_times-1]*output2.D[stimulus_times-1]/(output2.F[0]*output2.D[0]))
    #print(output2.D[stimulus_times])

    #print(EPSC_20hz)
    print(metrics.r2_score(EPSC_20hz, output1.EPSC))
    print(metrics.r2_score(EPSC_20hz, output2.EPSC))
    fig = plt.plot(range(n_pulses), output1.EPSC)
    plt.plot(range(n_pulses), output2.EPSC)
    plt.plot(range(n_pulses), EPSC_20hz, 'b.')
    plt.show()

    """TEST ONE METHOD AT TWO FREQUENCIES ACROSS A RANGE OF INPUTS"""
    # r = 10
    # stimulus_times1 = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int) #evenly spaced stimuli at r Hz starting at t = 1000/r ms
    # max_time = int(1000/r*(n_pulses+3))
    #
    # r = 20
    # stimulus_times2 = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int) #evenly spaced stimuli at r Hz starting at t = 1000/r ms
    # max_time = int(1000/r*(n_pulses+3))
    #
    # search k_o
    # output1 = method(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    # output2 = method(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    #
    """TEST ALL 4 METHODS ON THE SAME SET OF INPUTS"""
    # output1 = regular_train(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    # output2 = poisson_train(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    # output3 = DittmanRK1(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    # output4 = DittmanRK45(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    #
    # fig, axs = plt.subplots(4, sharex = True, sharey = True)
    # axs[0].plot(output1.times, -1*output1.EPSC_func/output1.EPSC_func[output1.stimulus_times[0]])
    # axs[1].plot(output2.times, -1*output2.EPSC_func/output2.EPSC_func[output2.stimulus_times[0]])
    # axs[2].plot(output3.times, -1*output3.EPSC_func/output3.EPSC_func[output3.stimulus_times[0]])
    # axs[3].plot(output3.times, -1*output4.EPSC_func/output4.EPSC_func[output4.stimulus_times[0]])
    #
    """TEST ONE METHOD ON 1 INPUT ACROSS A RANGE OF FREQUENCIES"""
    # EPSC_8_10_trend = []
    #
    # fig, axs = plt.subplots(4)
    #
    # hz_range = np.linspace(0.1,100,1000)
    # for i in hz_range:
    #
    #     stimulus_times = np.linspace(1000/i,(1000/i)*n_pulses,n_pulses, dtype = int)
    #     max_time = int(1000/i*(n_pulses+3))
    #
    #     output = method(stimulus_times, max_time, N_T, 0, 0.35, T_F, T_D, K_D, 0.7, 20, 0, 0.3, T_E)
    #
    #     EPSC_8_10_trend.append(sum(output.EPSC[-3:])/3)
    #
    #     if i == 1:
    #         axs[1].plot(output.times, -1*output.EPSC_func/output.EPSC_func[output.stimulus_times[0]])
    #         axs[1].set_xlim(0,max_time)
    #     elif i == 10:
    #         axs[2].plot(output.times, -1*output.EPSC_func/output.EPSC_func[output.stimulus_times[0]])
    #         axs[2].set_xlim(0,max_time)
    #     elif i == 100:
    #         axs[3].plot(output.times, -1*output.EPSC_func/output.EPSC_func[output.stimulus_times[0]])
    #         axs[3].set_xlim(0,max_time)
    #
    # axs[0].semilogx(hz_range, np.asarray(EPSC_8_10_trend)/.35)
    # axs[0].set_xlim(0.1,max(hz_range))

    """TEST ONE METHOD ON 3 INPUTS ACROSS A RANGE OF FREQUENCIES"""
    # analyticSol = []
    # analyticSol2 = []
    # analyticSol3 = []
    # EPSC_8_10_trend = []
    # EPSC_8_10_trend2 = []
    # EPSC_8_10_trend3 = []
    #
    # hz_range = np.geomspace(0.1,100,1000)
    # for i in hz_range:
    #
    #     stimulus_times = np.linspace(1000/i,1000/i*n_pulses,n_pulses, dtype=int)
    #     max_time = int(1000/i*(n_pulses+3))
    #
    #     output = method(stimulus_times, max_time, N_T, 0, 0.35, T_F, T_D, K_D, 0.7, 20, 0, 0.3, T_E)
    #     analyticSol.append(output.D_ss*output.F_ss/output.F_1)
    #     EPSC_8_10_trend.append((sum(output.EPSC[-3:])/3/output.EPSC[0]))
    #
    #     output2 = method(stimulus_times, max_time, N_T, 3.1, 0.05, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     analyticSol2.append(output2.D_ss*output2.F_ss/output2.F_1)
    #     EPSC_8_10_trend2.append((sum(output2.EPSC[-3:])/3/output2.EPSC[0]))
    #
    #     output3 = method(stimulus_times, max_time, N_T, 2.2, 0.24, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     analyticSol3.append(output3.D_ss*output3.F_ss/output3.F_1)
    #     EPSC_8_10_trend3.append((sum(output3.EPSC[-3:])/3/output3.EPSC[0]))
    #
    # fig = plt.semilogx(hz_range, analyticSol, color = 'red')
    # plt.semilogx(hz_range, analyticSol2, color = 'green')
    # plt.semilogx(hz_range, analyticSol3, color = 'blue')
    # plt.xlim(0.1,100)
    # plt.xlabel("Stimulus Rate (HZ)")
    # plt.ylabel("ESPC_8_10/EPSC_1")
    # plt.legend(('Climbing fiber', 'Parallel fiber', 'Schaffer collateral'))
    # plt.semilogx(hz_range, EPSC_8_10_trend, 'ro', markevery = 100)
    # plt.semilogx(hz_range, EPSC_8_10_trend2, 'go', markevery = 100)
    # plt.semilogx(hz_range, EPSC_8_10_trend3, 'bo', markevery = 100)

    """PLOT 3 SETS OF INPUTS WITH RK45"""
    # fig, axs = plt.subplots(3,3)
    #
    # output = DittmanRK45(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, 0.7, 20, 0, 0.3, T_E)
    #
    # output2 = DittmanRK45(stimulus_times, max_time, N_T, 3.1, 0.05, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    # output3 = DittmanRK45(stimulus_times, max_time, N_T, 2.2, 0.24, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    # axs[0, 0].plot(output.times, output.F)
    # axs[0, 0].set_title("Climbing fiber")
    # axs[0, 0].set_ylabel("F")
    # axs[0, 0].set_ylim(0,1)
    # axs[0, 0].set_xlim(0,max_time)
    #
    # axs[1,0].plot(output.times, output.D)
    # axs[1,0].set_ylabel("D")
    # axs[1,0].set_ylim(0,1)
    # axs[1,0].set_xlim(0,max_time)
    #
    # axs[2,0].plot(output.times, -1*output.EPSC_func/output.EPSC_func[output.stimulus_times[0]])
    # axs[2,0].set_xlabel('time (ms)')
    # axs[2,0].set_ylabel("Normalized ESPC")
    # #axs[2,0].set_ylim(-1,0)
    # axs[2,0].set_xlim(0,max_time)
    #
    # axs[0,1].plot(output.times, output2.F)
    # axs[0,1].set_title("Parallel fiber")
    # axs[0,1].set_ylabel("F2")
    # axs[0,1].set_ylim(0,1)
    # axs[0,1].set_xlim(0,max_time)
    #
    # axs[1,1].plot(output.times, output2.D)
    # axs[1,1].set_ylabel("D2")
    # axs[1,1].set_ylim(0,1)
    # axs[1,1].set_xlim(0,max_time)
    #
    # axs[2,1].plot(output.times, -1*output2.EPSC_func/output2.EPSC_func[output2.stimulus_times[0]])
    # axs[2,1].set_xlabel('time (ms)')
    # axs[2,1].set_ylabel("Normalized ESPC2")
    # #axs[2,1].set_ylim(-1,0)
    # axs[2,1].set_xlim(0,max_time)
    #
    # axs[0,2].plot(output.times, output3.F)
    # axs[0,2].set_title("Schaffer collateral")
    # axs[0,2].set_ylabel("F3")
    # axs[0,2].set_ylim(0,1)
    # axs[0,2].set_xlim(0,max_time)
    #
    # axs[1,2].plot(output.times, output3.D)
    # axs[1,2].set_ylabel("D3")
    # axs[1,2].set_ylim(0,1)
    # axs[1,2].set_xlim(0,max_time)
    #
    # axs[2,2].plot(output.times, -1*output3.EPSC_func/output3.EPSC_func[output3.stimulus_times[0]])
    # axs[2,2].set_xlabel('time (ms)')
    # axs[2,2].set_ylabel("Normalized ESPC3")
    # #axs[2,2].set_ylim(-1,0)
    # axs[2,2].set_xlim(0,max_time)
    #
    # plt.show()
