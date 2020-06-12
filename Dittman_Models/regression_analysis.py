import numpy as np
import pandas as pd
from DittmanModels import regular_train, poisson_train, DittmanRK1, DittmanRK45
from matplotlib import pyplot as plt

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
    F_1 = 0.3 #initial release probability
    T_F = 50 #decay constant for CaX_F
    T_D = 100 #decay constant for CaX_D
    K_D = 2 #affinity of CaX_D for site
    k_0 = 1 #initial recovery rate
    k_max = 20 #maximum recovery rate
    #K_F = 1 #affinity of CaX_F for site
    delta_F = 0 #amount by which CaX_F increases as a result of stimulus
    delta_D = 0 #amount by which CaX_D increases as a result of stimulus
    T_E = 2 #decay constant of simulated EPSCs

    #climbing fiber analysis
    #file in /home/evan/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx

    data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx', usecols = 'C,D', nrows = 20)
    data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 22, nrows = 20)
    data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 44, nrows = 20)
    data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 66, nrows = 20)

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

    SSregs = []
    best = [method, method, method, method, 0, 0]
    search_delta_D = np.linspace(0,1,1000)

    n_consider = 20
    for i in search_delta_D:
        output = method(stimulus_times1, max_time1, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, i, T_E)
        output2 = method(stimulus_times2, max_time2, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, i, T_E)
        output3 = method(stimulus_times3, max_time3, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, i, T_E)
        output4 = method(stimulus_times4, max_time4, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, i, T_E)

        SSreg = sum((output.EPSC[0:n_consider]/output.EPSC[0] - data_1hz.to_numpy()[0:n_consider,0])**2 + (output2.EPSC[0:n_consider]/output2.EPSC[0] - data_10hz.to_numpy()[0:n_consider,0])**2 + (output3.EPSC[0:n_consider]/output3.EPSC[0] - data_20hz.to_numpy()[0:n_consider,0])**2 + (output4.EPSC[0:n_consider]/output4.EPSC[0] - data_50hz.to_numpy()[0:n_consider,0])**2)/(4*n_consider) #average sum of squares deviation from the mean for each point

        if 1 - SSreg > best[-1]:
            best = [output, output2, output3, output4, i, 1 - SSreg]

        SSregs.append(1 - SSreg)


    print(best[-1])
    print(best[-2])

    #plot the best result
    fig, axs = plt.subplots(4)
    axs[0].plot(best[0].times, -1*best[0].EPSC_func/best[0].EPSC_func[best[0].stimulus_times[0]])
    axs[0].errorbar(best[0].stimulus_times, -1*data_1hz.to_numpy()[:,0], yerr = data_1hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)

    axs[1].plot(best[1].times, -1*best[1].EPSC_func/best[1].EPSC_func[best[1].stimulus_times[0]])
    axs[1].errorbar(best[1].stimulus_times, -1*data_10hz.to_numpy()[:,0], yerr = data_10hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)

    axs[2].plot(best[2].times, -1*best[2].EPSC_func/best[2].EPSC_func[best[2].stimulus_times[0]])
    axs[2].errorbar(best[2].stimulus_times, -1*data_20hz.to_numpy()[:,0], yerr = data_20hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)

    axs[3].plot(best[3].times, -1*best[3].EPSC_func/best[3].EPSC_func[best[3].stimulus_times[0]])
    axs[3].errorbar(best[3].stimulus_times, -1*data_50hz.to_numpy()[:,0], yerr = data_50hz.to_numpy()[:,1], fmt = '.r', ecolor = 'black', elinewidth = 10)

#testing outputs

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
