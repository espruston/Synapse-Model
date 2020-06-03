import matplotlib.pyplot as plt
import numpy as np
import sys
from math import e

def calculate_CaX_F(T_F, delta_F, stimuli):

    CaX_F = np.zeros(max_time+1) #preallocate vector of CaX_F values from t=0->max_time
    CaX_F = CaX_F + 1e-30
    for t in times:
        dCaX_F_dt = -1*CaX_F[t-1]/T_F + delta_F*stimuli[t-1] #exponential decay with characteristic time T_F after a jump of delta_F after an action potential (Eq. 3)
        CaX_F[t] = CaX_F[t-1] + dCaX_F_dt #CaX_F[0] should be 0 so that no facilitation occurs before an action potential

    return CaX_F

def calculate_F(CaX_F, K_F, F_1, times):

    F = np.zeros(max_time+1) #preallocate vector of F values from t=0->max_time
    F[0] = F_1 #F = F_1 until a stimulus occurs
    for t in times:
        F[t] = F_1 + (1 - F_1)/(1 + K_F/CaX_F[t]) #(Eq. 2a)

    return F

def calculate_CaX_D(T_D, delta_D, stimuli):

    CaX_D = np.zeros(max_time+1) #preallocate vector of CaX_D values from t=0->max_time
    CaX_D = CaX_D + 1e-19
    for t in times:
        dCaX_D_dt = -1*CaX_D[t-1]/T_D + delta_D*stimuli[t-1] #exponential decay with characteristic time T_D after a jump of delta_D after an action potential (Eq. 12)
        CaX_D[t] = CaX_D[t-1] + dCaX_D_dt #CaX_D[0] should be 0 so that np Ca dependent recovery occurs before an action potential

    return CaX_D

def calculate_k_recov(CaX_D, K_D, k_0, k_max):

    k_recov = np.zeros(max_time+1) #preallocate values of k_recov
    k_recov[0] = k_0
    for t in times:
        k_recov[t] = (k_max - k_0)/(1 + K_D/CaX_D[t]) + k_0

    return k_recov

def calculate_D(CaX_D, k_recov, F, stimuli):

    D = np.ones(max_time+1) #preallocate vector of D values from t=0->max_time (D=1 @ t=0)
    for t in times:
        dD_dt = (1 - D[t-1])*k_recov[t-1] - D[t-1]*F[t-1]*stimuli[t-1] #scale k_recov to msec
        #dD_dt = (1 - D[t-1])*k_recov/1000 - D[t-1]*F[t-1]*stimuli[t-1]
        D[t] = D[t-1] + dD_dt

    return D

def calculate_EPSC(N_T, F, D, T_E):

    alpha_1 = (times*e/T_E)*e**(-1*times/T_E) #reference alpha function
    alpha = np.zeros(max_time) #functional alpha function

    for stimulus in stimulus_times:
        alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (F[stimulus-1] * D[stimulus-1] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

    EPSC =  alpha * N_T

    return EPSC, alpha

if __name__ == "__main__": #when called by python interpreter as "python DittmanModel.py"

    # if len(sys.argv) != 11: #incorrect number of arguments leads to exit and error message
    #     print 'ERROR! Call program as: DittmanModel.py F_1 T_D K_D N_T K_F T_F k_0 k_max delta_F delta_D'
    #     sys.exit(2)
    #
    # else: #assign passed values to parameters
    #     F_1 = int(sys.argv[1]) #initial probability of release (0->1)
    #     T_D = int(sys.argv[2]) #decay time constant of CaX_D after an action potential (units of sec)
    #     K_D = int(sys.argv[3]) #affinity of CaX_D for release site
    #     N_T = int(sys.argv[4]) #total number of release sites
    #     K_F = int(sys.argv[5]) #affinity of CaX_F for release site
    #     T_F = int(sys.argv[6]) #decay time constant of CaX_F after an action potential
    #     k_0 = int(sys.argv[7]) #baseline rate of recovery from refractory state (units of sec^-1)
    #     k_max = int(sys.argv[8]) #maximal rate of recovery from refractory state (units of sec^-1)
    #     delta_F = int(sys.argv[9]) #amount by which F increases as a result of an action potential
    #     delta_D = int(sys.argv[10]) #amount by which D increases as a result of an action potential

        args = [0.35, 50, 2, 1, 2, 1e30, 0.7, 20, 0, 1, 2]
        F_1 = args[0]
        T_D = args[1]
        K_D = args[2]
        N_T = args[3]
        K_F = args[4]
        T_F = args[5]
        k_0 = args[6]/1000 #timesteps are msec
        k_max  = args[7]/1000 #timestamps are msec
        delta_F = args[8]
        delta_D = args[9]
        T_E = args[10]

        max_time = 3000 #3000 msec
        times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        stimulus_times = [20,40,60,80,100,120,140,160,180,200,220] #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector

        stimuli = np.zeros(max_time)
        for stimulus in stimulus_times:
            stimuli[stimulus] = 1 #stiumli vector now is 1 for times when a stimulus occurs and 0 if it doesn't

        CaX_F = calculate_CaX_F(T_F, delta_F, stimuli)
        F = calculate_F(CaX_F, K_F, F_1, times)
        CaX_D = calculate_CaX_D(T_D, delta_D, stimuli)
        k_recov = calculate_k_recov(CaX_D, K_D, k_0, k_max)
        D = calculate_D(CaX_D, k_recov, F, stimuli)
        EPSC, alpha = calculate_EPSC(N_T, F, D, T_E)

        fig, axs = plt.subplots(3, sharex='col')

        axs[0].plot(times,-1*EPSC)
        axs[1].plot(times,F[1:3001])
        axs[2].plot(times,D[1:3001])
        plt.show()
