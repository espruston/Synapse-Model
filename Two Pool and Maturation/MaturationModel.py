import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from math import e

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

    n_sites = 1
    n_pulses = 20
    r = 1
    delta_t = 1/r

    #release probabilities
    #p_immature = 0.1
    #p_mature = .5
    p_facilitated = .9

    #time constants
    T_refill = 1e-30
    #T_maturation = 1e30
    T_facilitation = 1e30

    #pool sizes
    n_mature = [n_sites]
    n_immature = [0]
    n_facilitated = [0]
    n_empty = [0]

    #release list preallocation
    n_released_facilitated = []
    n_released_mature = []
    n_released_immature = []
    n_released_tot = []

    EPSC = []
    stdevs = []

    best = [1e30]
    search_T_maturation = np.linspace(0.01,10.01,110)
    search_p_mature = np.linspace(0.02,0.9,89)
    search_p_immature = np.linspace(0,0.22,100)

    n_consider = 20

    for T_maturation in search_T_maturation:
        for p_mature in search_p_mature:
            for p_immature in search_p_immature:

                n_mature = [n_sites]
                n_immature = [0]
                n_facilitated = [0]
                n_empty = [0]

                n_released_facilitated = []
                n_released_mature = []
                n_released_immature = []
                n_released_tot = []

                EPSC = []

                for i in range(n_pulses):
                    #simulate a pulse, then pass through delta_t before simulating next pulse
                    n_released_facilitated.append(n_facilitated[i]*p_facilitated)
                    n_released_mature.append(n_mature[i]*p_mature)
                    n_released_immature.append(n_immature[i]*p_immature)

                    n_released_tot.append(n_released_facilitated[i] + n_released_mature[i] + n_released_immature[i])
                    EPSC.append(n_released_tot[i])

                    n_empty.append((n_empty[i] + n_released_tot[i])*(e**(-delta_t/T_refill))) #after vesicles are released, there are n_released_tot additional empty sites which refill with characteristic time T_refill over delta_t
                    n_mature.append(n_mature[i] - n_released_mature[i] + (n_immature[i] - n_released_immature[i])*(1-e**(-delta_t/T_maturation))) #some of the mature pool is released, some is replenished by the immature pool that remains after the pulse and some is replenished by the new members joining the immature pool first
                    #what equation defines the third portion of replenishment
                    n_immature.append((n_immature[i] - n_released_immature[i])*(e**(-delta_t/T_maturation)) + (n_empty[i] + n_released_tot[i])*(1-e**(-delta_t/T_refill)))  #some immature vesicles mature and leave the immature pool, some vesicles join the immature pool over delta_t
                    n_facilitated.append(0) #no facilitation is occuring

                stdev = np.sqrt(sum((np.asarray(EPSC[0:n_consider])/EPSC[0] - EPSC_1hz[0:n_consider])**2)/n_consider)

                if stdev < best[-1]:
                    best = [EPSC, T_maturation, p_mature, p_immature, stdev]

                stdevs.append(stdev)

    stimulus_times = np.linspace(1000/r,(1000/r)*n_pulses,n_pulses, dtype = int)
    max_time = int(1000/r*(n_pulses+3))
    times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

    alpha_1 = (times*e/2)*e**(-1*times/2) #reference alpha function
    alpha = np.zeros(max_time) #functional alpha function

    for i in range(n_pulses):
        stimulus = stimulus_times[i]
        alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (best[0][i] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

    fig = plt.plot(times, -1*alpha/alpha[stimulus_times[0]])
    plt.errorbar(stimulus_times, -1*EPSC_1hz, yerr = stdev_1hz, fmt = '.r', ecolor = 'black', elinewidth = 10 )
    plt.show()
