import numpy as np
import pandas as pd
import os
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib import pyplot as plt
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
    delta_D = 0 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    #climbing fiber analysis
    #file in /home/evan/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx

    data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', nrows = 20)
    data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 22, nrows = 20)
    data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 44, nrows = 20)
    data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 66, nrows = 20)

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

    # best = [method, method, method, method, 0, 1e30]
    best = [0]
    # search_K_D = np.linspace(0,100,11)
    # search_T_D = np.linspace(1,100,100)
    # search_k_max = np.linspace(k_0,99+k_0,100)

    stdevs = []
    r_squareds = []
    #ss_weight = .5

    EPSCs1 = []
    EPSCs10 = []
    EPSCs20 = []
    EPSCs50 = []
    EPSCss = []

    n_consider = 20

    data_1hz_np = data_1hz.to_numpy()[0:n_consider,0]
    data_10hz_np = data_10hz.to_numpy()[0:n_consider,0]
    data_20hz_np = data_20hz.to_numpy()[0:n_consider,0]
    data_50hz_np = data_50hz.to_numpy()[0:n_consider,0]

    stdev_1hz = data_1hz.to_numpy()[:,1]
    stdev_10hz = data_10hz.to_numpy()[:,1]
    stdev_20hz = data_20hz.to_numpy()[:,1]
    stdev_50hz = data_50hz.to_numpy()[:,1]

    ss_stdev_1hz = data_1hz.to_numpy()[n_pulses-1,1]*5
    ss_stdev_10hz = data_10hz.to_numpy()[n_pulses-1,1]*5
    ss_stdev_20hz = data_20hz.to_numpy()[n_pulses-1,1]*5
    ss_stdev_50hz = data_50hz.to_numpy()[n_pulses-1,1]*5

    ss_1hz = data_1hz.to_numpy()[n_pulses-1, 0]
    ss_10hz = data_10hz.to_numpy()[n_pulses-1, 0]
    ss_20hz = data_20hz.to_numpy()[n_pulses-1, 0]
    ss_50hz = data_50hz.to_numpy()[n_pulses-1, 0]

    ss_range_1hz = [ss_1hz - ss_stdev_1hz, ss_1hz + ss_stdev_1hz]
    ss_range_10hz = [ss_10hz - ss_stdev_10hz, ss_10hz + ss_stdev_10hz]
    ss_range_20hz = [ss_20hz - ss_stdev_20hz, ss_20hz + ss_stdev_20hz]
    ss_range_50hz = [ss_50hz - ss_stdev_50hz, ss_50hz + ss_stdev_50hz]

    # #TEST AT 10 AND 20 hz, FIT k_0
    # search_k_0 = np.linspace(0.01,20.01,200)
    # for k_0 in search_k_0:
    #     output1 = method(stimulus_times2, max_time2, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     output2 = method(stimulus_times3, max_time3, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    #     r_squared1 = 1 - (sum((output1.EPSC[0:n_consider]/output1.EPSC[0] - data_10hz_np)**2)/sum((data_10hz_np - np.average(data_10hz_np))**2))
    #
    #     r_squared2 = 1 - (sum((output2.EPSC[0:n_consider]/output2.EPSC[0] - data_20hz_np)**2)/sum((data_20hz_np - np.average(data_20hz_np))**2))
    #
    #     r_squareds.append([r_squared1, r_squared2, (r_squared1 + r_squared2)/2])
    #
    #     EPSCs10.append(output1.EPSC[0:n_consider]/output1.EPSC[0])
    #     EPSCs20.append(output2.EPSC[0:n_consider]/output2.EPSC[0])
    #
    #     if (r_squared1 + r_squared2)/2 > best[-1]:
    #         best = [k_0, EPSCs10[-1], EPSCs20[-1], r_squared1, r_squared2, (r_squared1 + r_squared2)/2]
    #
    # pathset = os.path.expanduser(r"~/Dropbox/Work/Jackman Lab/Modeling/200701_1_10hz_20hz_0F0D")
    # np.savez(pathset, search_k_0, best[0], r_squareds, best[-3:], best[1:2], EPSCs10, EPSCs20)
    #
    # print("best average R squared =", best[-1]) #avg r_squared
    #
    # fig, axs = plt.subplots(2)
    #
    # axs[0].errorbar(stimulus_times1-stimulus_times1[0], data_10hz_np, yerr = stdev_10hz, label = "data")
    # axs[0].plot(stimulus_times1-stimulus_times1[0], best[1], label = "sim")
    # axs[0].set_xlim(0,max_time1)
    # axs[0].set_ylim(0,1.25)
    # axs[0].set_xlabel("t (ms)")
    # axs[0].set_ylabel("$EPSC(t)/EPSC_{1}$")
    # axs[0].set_title("10 Hz Stim, $R^{2}=$" + str(best[-3]) + " for $k_{0}=$" + str(best[0]))
    # axs[0].legend()
    #
    # axs[1].errorbar(stimulus_times2-stimulus_times2[0], data_20hz_np, yerr = stdev_20hz, label = "data")
    # axs[1].plot(stimulus_times2-stimulus_times2[0], best[2], label = "sim")
    # axs[1].set_xlim(0,max_time2)
    # axs[1].set_ylim(0,1.25)
    # axs[1].set_xlabel("t (ms)")
    # axs[1].set_ylabel("$EPSC(t)/EPSC_{1}$")
    # axs[1].set_title("20 Hz Stim, $R^{2}=$" + str(best[-2]) + " for $k_{0}=$" + str(best[0]))
    # axs[1].legend()

    #TEST AT 1, 10, 20 AND 50 hz, FIT k_0 with 0 facil and Ca dependent recovery
    # search_k_0 = np.linspace(0.01,20.01,200)
    # for k_0 in search_k_0:
    #     output1 = method(stimulus_times1, max_time1, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     output2 = method(stimulus_times2, max_time2, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     output3 = method(stimulus_times3, max_time3, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #     output4 = method(stimulus_times4, max_time4, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E)
    #
    #     r_squared1 = 1 - (sum((output1.EPSC[0:n_consider]/output1.EPSC[0] - data_1hz_np)**2)/sum((data_1hz_np - np.average(data_1hz_np))**2))
    #
    #     r_squared2 = 1 - (sum((output2.EPSC[0:n_consider]/output2.EPSC[0] - data_10hz_np)**2)/sum((data_10hz_np - np.average(data_10hz_np))**2))
    #
    #     r_squared3 = 1 - (sum((output3.EPSC[0:n_consider]/output3.EPSC[0] - data_20hz_np)**2)/sum((data_20hz_np - np.average(data_20hz_np))**2))
    #
    #     r_squared4 = 1 - (sum((output4.EPSC[0:n_consider]/output4.EPSC[0] - data_50hz_np)**2)/sum((data_50hz_np - np.average(data_50hz_np))**2))
    #
    #     r_squareds.append([r_squared1, r_squared2, r_squared3, r_squared4, (r_squared1 + r_squared2 + r_squared3 + r_squared4)/4])
    #
    #     EPSCs1.append(output1.EPSC[0:n_consider]/output1.EPSC[0])
    #     EPSCs10.append(output2.EPSC[0:n_consider]/output2.EPSC[0])
    #     EPSCs20.append(output3.EPSC[0:n_consider]/output3.EPSC[0])
    #     EPSCs50.append(output4.EPSC[0:n_consider]/output4.EPSC[0])
    #
    #     if (r_squared1 + r_squared2 + r_squared3 + r_squared4)/4 > best[-1]:
    #         best = [k_0, EPSCs1[-1], EPSCs10[-1], EPSCs20[-1], EPSCs50[-1], r_squared1, r_squared2, r_squared3, r_squared4, (r_squared1 + r_squared2 + r_squared3 + r_squared4)/4]
    #
    # pathset = os.path.expanduser(r"~/Dropbox/Work/Jackman Lab/Modeling/200701_1_1hz_10hz_20hz_50hz_0F0D")
    # np.savez(pathset, best[0], best[1:4], best[-5:], search_k_0,  r_squareds, EPSCs1, EPSCs10, EPSCs20, EPSCs50)
    #
    # print("best average R squared =", best[-1]) #avg r_squared
    #
    # fig, axs = plt.subplots(2,2)
    #
    # axs[0,0].errorbar(stimulus_times1-stimulus_times1[0], data_1hz_np, yerr = stdev_1hz, label = "data", fmt = '.r')
    # axs[0,0].plot(stimulus_times1-stimulus_times1[0], best[1], label = "sim")
    # axs[0,0].set_xlim(0,max_time1)
    # axs[0,0].set_ylim(0,1.25)
    # axs[0,0].set_xlabel("t (ms)")
    # axs[0,0].set_ylabel("$EPSC(t)/EPSC_{1}$")
    # axs[0,0].set_title("1 Hz Stim, $R^{2}=$" + str(best[-5]) + " for $k_{0}=$" + str(best[0]))
    # axs[0,0].legend()
    #
    # axs[1,0].errorbar(stimulus_times2-stimulus_times2[0], data_10hz_np, yerr = stdev_10hz, label = "data", fmt = '.r')
    # axs[1,0].plot(stimulus_times2-stimulus_times2[0], best[2], label = "sim")
    # axs[1,0].set_xlim(0,max_time2)
    # axs[1,0].set_ylim(0,1.25)
    # axs[1,0].set_xlabel("t (ms)")
    # axs[1,0].set_ylabel("$EPSC(t)/EPSC_{1}$")
    # axs[1,0].set_title("10 Hz Stim, $R^{2}=$" + str(best[-4]) + " for $k_{0}=$" + str(best[0]))
    # axs[1,0].legend()
    #
    # axs[0,1].errorbar(stimulus_times3-stimulus_times3[0], data_20hz_np, yerr = stdev_20hz, label = "data", fmt = '.r')
    # axs[0,1].plot(stimulus_times3-stimulus_times3[0], best[3], label = "sim")
    # axs[0,1].set_xlim(0,max_time3)
    # axs[0,1].set_ylim(0,1.25)
    # #axs[2].set_xlabel("t (ms)")
    # axs[0,1].set_ylabel("$EPSC(t)/EPSC_{1}$")
    # axs[0,1].set_title("20 Hz Stim, $R^{2}=$" + str(best[-3]) + " for $k_{0}=$" + str(best[0]))
    # axs[0,1].legend()
    #
    # axs[1,1].errorbar(stimulus_times4-stimulus_times4[0], data_50hz_np, yerr = stdev_50hz, label = "data", fmt = '.r')
    # axs[1,1].plot(stimulus_times4-stimulus_times4[0], best[4], label = "sim")
    # axs[1,1].set_xlim(0,max_time4)
    # axs[1,1].set_ylim(0,1.25)
    # axs[1,1].set_xlabel("t (ms)")
    # axs[1,1].set_ylabel("$EPSC(t)/EPSC_{1}$")
    # axs[1,1].set_title("50 Hz Stim, $R^{2}=$" + str(best[-2]) + " for $k_{0}=$" + str(best[0]))
    # axs[1,1].legend()

    # #TEST WITHOUT FACIL OR RECOVERY MECHANISM
    # for T_0 in np.linspace(0.1,10.1, 100):
    #     for F in np.linspace(0.5,1, 60):
    #         filled_pre = [1]
    #         filled_post = []
    #         delta_t = 20
    #         for j in range(n_consider):
    #             filled_post.append(filled_pre[j]*(1-F))
    #             filled_pre.append(filled_post[-1] + (1-filled_post[-1])*e**(-delta_t/T_0))
    #         stdev = np.sqrt(sum((np.asarray(filled_pre[0:n_consider]) - data_50hz_np)**2))
    #         if stdev < best[-1]:
    #             best = [np.asarray(filled_pre), T_0, F, stdev]
    #
    # fig = plt.scatter(range(n_consider), best[0][0:n_consider]*-1)
    # plt.errorbar(range(n_consider), -1*data_50hz_np, yerr = stdev_50hz[0:n_consider], fmt = '.r', ecolor = 'black', elinewidth = 10)
    # plt.show()

    #SEARCH K_D T_D and k_max
    # for i in search_K_D:
    #     for j in search_T_D:
    #         for k in search_k_max:
    #             output1 = method(stimulus_times1, max_time1, N_T, roh, F_1, T_F, j, i, k_0, k, delta_F, delta_D, T_E)
    #             output2 = method(stimulus_times2, max_time2, N_T, roh, F_1, T_F, j, i, k_0, k, delta_F, delta_D, T_E)
    #             output3 = method(stimulus_times3, max_time3, N_T, roh, F_1, T_F, j, i, k_0, k, delta_F, delta_D, T_E)
    #             output4 = method(stimulus_times4, max_time4, N_T, roh, F_1, T_F, j, i, k_0, k, delta_F, delta_D, T_E)
    #
    #             stdev = np.sqrt(sum((output4.EPSC[0:n_consider]/output4.EPSC[0] - data_50hz_np)**2)/n_consider)
    #
    #             # stdev = np.sqrt(sum((output1.EPSC[0:n_consider]/output1.EPSC[0] - data_1hz_np)**2 + (output2.EPSC[0:n_consider]/output2.EPSC[0] - data_10hz_np)**2 + (output3.EPSC[0:n_consider]/output3.EPSC[0] - data_20hz_np)**2 + (output4.EPSC[0:n_consider]/output4.EPSC[0] - data_50hz_np)**2)/(4*n_consider)) #average sqrt of sum of squares deviation from the mean for each point
    #
    #             EPSCs1.append(output1.EPSC[0:n_consider])
    #             EPSCs10.append(output2.EPSC[0:n_consider])
    #             EPSCs20.append(output3.EPSC[0:n_consider])
    #             EPSCs50.append(output4.EPSC[0:n_consider])
    #
    #             # ssdev = np.sqrt(((output1.EPSC[n_pulses-1]/output1.EPSC[0] - ss_1hz)**2 + (output2.EPSC[n_pulses-1]/output2.EPSC[0] - ss_10hz)**2 + (output3.EPSC[n_pulses-1]/output3.EPSC[0] - ss_20hz)**2 + (output4.EPSC[n_pulses-1]/output4.EPSC[0] - ss_50hz)**2)/4)
    #             # EPSCss.append(output1.EPSC[n_pulses-1])
    #             # EPSCss.append(output2.EPSC[n_pulses-1])
    #             # EPSCss.append(output3.EPSC[n_pulses-1])
    #             # EPSCss.append(output4.EPSC[n_pulses-1])
    #
    #
    #             # if (ss_range_1hz[0] <= output1.EPSC_norm_ss <= ss_range_1hz[1]) and (ss_range_10hz[0] <= output2.EPSC_norm_ss <= ss_range_10hz[1]) and (ss_range_20hz[0] <= output3.EPSC_norm_ss <= ss_range_20hz[1]) and (ss_range_50hz[0] <= output4.EPSC_norm_ss <= ss_range_50hz[1]) and (stdev < best[-1]): #all steady state values are within 5 SEM
    #             if stdev < best[-1]:
    #                 best = [output1, output2, output3, output4, i, j, k, stdev]
    #
    #             stdevs.append(stdev)

    # print(best[-1])
    # print(best[-2])
    # print(best[-3])
#    print(best[-4])

    # EPSCs1 = np.asarray(EPSCs1).reshape(n_consider, len(search_K_D), len(search_T_D), len(search_k_max))
    # EPSCs10 = np.asarray(EPSCs10).reshape(n_consider, len(search_K_D), len(search_T_D), len(search_k_max))
    # EPSCs20 = np.asarray(EPSCs20).reshape(n_consider, len(search_K_D), len(search_T_D), len(search_k_max))
    # EPSCs50 = np.asarray(EPSCs50).reshape(n_consider, len(search_K_D), len(search_T_D), len(search_k_max))
    # stdevs = np.asarray(stdevs).reshape(len(search_K_D), len(search_T_D), len(search_k_max))

#     np.savez('/home/evan/Documents/Work/JackmanLab/EPSCs.200618_1', np.asarray([best[-4], best[-3], best[-2]]), search_K_D, search_T_D, search_k_max, stdevs, EPSCs1, EPSCs10, EPSCs20, EPSCs50)
#     #plot best 50hz fit
#     # fig = plt.plot(np.linspace(1, n_consider, n_consider), best[3].EPSC[0:n_consider]/best[3].EPSC[0])
#     # plt.errorbar(np.linspace(1, n_consider, n_consider), data_50hz.to_numpy()[0:n_consider,0], yerr = data_50hz.to_numpy()[0:n_consider,1], fmt = '.r', ecolor = 'black', elinewidth = 10)
#
#     #plot the best result
#     fig, axs = plt.subplots(4)
#     axs[0].plot(best[0].times, -1*best[0].EPSC_func/best[0].EPSC_func[best[0].stimulus_times[0]])
#     axs[0].errorbar(best[0].stimulus_times, -1*data_1hz.to_numpy()[:,0], yerr = stdev_1hz, fmt = '.r', ecolor = 'black', elinewidth = 10)
#
#     axs[1].plot(best[1].times, -1*best[1].EPSC_func/best[1].EPSC_func[best[1].stimulus_times[0]])
#     axs[1].errorbar(best[1].stimulus_times, -1*data_10hz.to_numpy()[:,0], yerr = stdev_10hz, fmt = '.r', ecolor = 'black', elinewidth = 10)
#
#     axs[2].plot(best[2].times, -1*best[2].EPSC_func/best[2].EPSC_func[best[2].stimulus_times[0]])
#     axs[2].errorbar(best[2].stimulus_times, -1*data_20hz.to_numpy()[:,0], yerr = stdev_20hz, fmt = '.r', ecolor = 'black', elinewidth = 10)
#
#     axs[3].plot(best[3].times, -1*best[3].EPSC_func/best[3].EPSC_func[best[3].stimulus_times[0]])
#     axs[3].errorbar(best[3].stimulus_times, -1*data_50hz.to_numpy()[:,0], yerr = stdev_50hz, fmt = '.r', ecolor = 'black', elinewidth = 10)
# #testing outputs

    # print(data_1hz)
    # print(data_10hz)
    # print(data_20hz)
    # print(data_50hz)

    # print(np.asarray(output.EPSC[-1])/output.EPSC[0])
    # print(np.asarray(output2.EPSC[-1])/output2.EPSC[0])
    # print(np.asarray(output3.EPSC[-1])/output3.EPSC[0])
    # print(np.asarray(output4.EPSC[-1])/output4.EPSC[0])

    # fig, axs = plt.subplots(4)
    # axs[0].plot(output.times, -1*output.EPSC_func/output.EPSC_func[output.stimulus_times[0]])
    # axs[0].errorbar(output.stimulus_times, -1*data_1hz.to_numpy()[:,0], yerr = data_1hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)
    #
    # axs[1].plot(output2.times, -1*output2.EPSC_func/output2.EPSC_func[output2.stimulus_times[0]])
    # axs[1].errorbar(output2.stimulus_times, -1*data_10hz.to_numpy()[:,0], yerr = data_10hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)
    #
    # axs[2].plot(output3.times, -1*output3.EPSC_func/output3.EPSC_func[output3.stimulus_times[0]])
    # axs[2].errorbar(output3.stimulus_times, -1*data_20hz.to_numpy()[:,0], yerr = data_20hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)
    #
    # axs[3].plot(output4.times, -1*output4.EPSC_func/output4.EPSC_func[output4.stimulus_times[0]])
    # axs[3].errorbar(output4.stimulus_times, -1*data_50hz.to_numpy()[:,0], yerr = data_50hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)
    plt.show()
