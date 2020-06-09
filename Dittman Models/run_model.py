import numpy as np
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib import pyplot as plt

if __name__ == "__main__":

    METHODS = {'regular': regular_train,
               'poisson': poisson_train,
               'RK1': DittmanRK1,
               'RK45': DittmanRK45}

    method = 'RK45'
    method = METHODS[method]

    #arguments passed to simulator
    n_pulses = 10 #only needed for regular train
    r = 50 #only needed for regular train
    if method == METHODS['regular']:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int) #evenly spaced stimuli at r Hz starting at t = 1000/r ms
    else:
        stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
        #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector [0,1,2,...,max_time]
    max_time = 320 #used in plotting
    N_T = 1 #number of total release sites
    roh = 2.2 #EPSC2/EPSC1
    F_1 = 0.24 #initial facilitation
    T_F = 100 #decay constant for CaX_F
    T_D = 50 #decay constant for CaX_D
    K_D = 2 #affinity of CaX_D for site
    k_0 = 2 #initial recovery rate
    k_max = 30 #maximum recovery rate
    #K_F = 1 #affinity of CaX_F for site
    delta_F = 1 #amount by which CaX_F increases as a result of stimulus
    delta_D = 1 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    output1 = regular_train(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    output2 = poisson_train(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    output3 = DittmanRK1(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)

    output4 = DittmanRK45(stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)


    fig, axs = plt.subplots(4, sharex = True, sharey = True)
    axs[0].plot(output1.times, -1*output1.EPSC_func/output1.EPSC_func[output1.stimulus_times[0]])
    axs[1].plot(output2.times, -1*output2.EPSC_func/output2.EPSC_func[output2.stimulus_times[0]])
    axs[2].plot(output3.times, -1*output3.EPSC_func/output3.EPSC_func[output3.stimulus_times[0]])
    axs[3].plot(output3.times, -1*output4.EPSC_func/output4.EPSC_func[output4.stimulus_times[0]])
    # EPSC_8_10_trend = []
    # EPSC_8_10_trend_2 = []
    # EPSC_8_10_trend_3 = []
    #
    # hz_range = np.linspace(1,100,100,dtype=int)
    # for i in hz_range:
    #     print(i)
    #
    #     stimulus_times = np.linspace(1000/i,(1000/i)*n_pulses,n_pulses, dtype = int)
    #     max_time = int(1000/i*n_pulses)
    #
    #     output = method(stimulus_times, max_time, N_T, 0, 0.35, T_F, T_D, K_D, k_0, k_max, 0, delta_D, T_E)
    #
    #     EPSC_8_10_trend.append(sum(output.EPSC_func[stimulus_times[-3:]-1])/3)
    #
    #     output2 = method(stimulus_times, max_time, N_T, 3.1, 0.05, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     EPSC_8_10_trend_2.append(sum(output2.EPSC_func[stimulus_times[-3:]-1])/3)
    #
    #     output3 = method(stimulus_times, max_time, N_T, 2.2, 0.24, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     EPSC_8_10_trend_3.append(sum(output3.EPSC_func[stimulus_times[-3:]-1])/3)
    #
    # fig, axs = plt.subplots(3)
    # axs[0].plot(hz_range, np.asarray(EPSC_8_10_trend)/.35)
    # axs[1].plot(hz_range, np.asarray(EPSC_8_10_trend_2)/.05)
    # axs[2].plot(hz_range, np.asarray(EPSC_8_10_trend_3)/.24)

    # fig, axs = plt.subplots(3,3)
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
    # axs[2,0].plot(output.times, -1*output.EPSC_func/max(output.EPSC_func))
    # axs[2,0].set_xlabel('time (ms)')
    # axs[2,0].set_ylabel("Normalized ESPC")
    # axs[2,0].set_ylim(-1,0)
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
    # axs[2,1].plot(output.times, -1*output2.EPSC_func/max(output2.EPSC_func))
    # axs[2,1].set_xlabel('time (ms)')
    # axs[2,1].set_ylabel("Normalized ESPC2")
    # axs[2,1].set_ylim(-1,0)
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
    # axs[2,2].plot(output.times, -1*output3.EPSC_func/max(output3.EPSC_func))
    # axs[2,2].set_xlabel('time (ms)')
    # axs[2,2].set_ylabel("Normalized ESPC3")
    # axs[2,2].set_ylim(-1,0)
    # axs[2,2].set_xlim(0,max_time)

    plt.show()
