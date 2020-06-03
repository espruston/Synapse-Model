import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from math import e

def dCaX_F_dt(t, CaX_F):

    return -1*CaX_F/T_F + delta_F*(int(t) in stimulus_times) #exponential decay with characteristic time T_F after a jump of delta_F after an action potential (Eq. 3)

def dCaX_D_dt(t, CaX_D):

    return -1*CaX_D/T_D + delta_D*(int(t) in stimulus_times) #exponential decay with characteristic time T_D after a jump of delta_D after an action potential (Eq. 12)

def dD_dt(t, D):

    CaX_F_temp = CaX_F.sol(t) + 1e-30*(CaX_F.sol(t) == 0) #sets F to F_1 when CaX_F(t) = 0
    CaX_D_temp = CaX_D.sol(t) + 1e-30*(CaX_D.sol(t) == 0) #sets k_recov to k_0 when CaX_D(t) = 0

    k_recov_temp = (k_max - k_0)/(1 + K_D/CaX_D_temp) + k_0
    k_recov.append(k_recov_temp)

    return (1 - D)*k_recov_temp - D*F[0][int(t)-1]*(int(t) in stimulus_times) #eq 13

if __name__ == "__main__": #when called by python interpreter as "python DittmanModelRK45.py"

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

        args = [0.35, 50, 2, 1, 2, 100, 0.7, 20, 1, 1, 2]
        F_1 = args[0]
        T_D = args[1]
        K_D = args[2]
        N_T = args[3]
        K_F = args[4]
        T_F = args[5]
        k_0 = args[6]/1000 #timesteps are msec
        k_max  = args[7]/1000 #timesteps are msec
        delta_F = args[8]
        delta_D = args[9]
        T_E = args[10]

        max_time = 1300 #number in msec
        times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        stimulus_times = [10,30,50,70,90,110,130,150,170,190,1010,1030,1050,1070,1090,1110,1130,1150,1170,1190] #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector
        max_step_size = 1

        CaX_F = solve_ivp(dCaX_F_dt, [0, max_time], [0], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = max_step_size) #solve dCaX_F_dt using RK45 starting from first stimulus w max step size 1

        CaX_D = solve_ivp(dCaX_D_dt, [0, max_time], [0], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = max_step_size) #solve dCaX_D_dt using RK45 starting from first stimulus w max step size 1

        F = F_1 + (1 - F_1)/(1 + K_F/(CaX_F.sol(times) + 1e-30*(CaX_F.sol(times) == 0)))
        k_recov = []

        D = solve_ivp(dD_dt, [0, max_time], [1], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = max_step_size) #solve dD_dt using RK45 starting from first stimulus w max step size 1

        alpha_1 = (times*e/T_E)*e**(-1*times/T_E) #reference alpha function
        alpha = np.zeros(max_time) #functional alpha function is preallocated as all zeros from t=1->max_time

        for stimulus in stimulus_times: #for each stimulus, scale it by F*D, shift it to the stimulus time, and add it to the functional alpha function
            alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (F[0][stimulus-1] * D.sol(stimulus) * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        EPSC =  -1*alpha * N_T

        #plot variables as a function of time

        fig, axs = plt.subplots(3, sharex='col')

        axs[0].plot(times, F[0])
        axs[0].set_ylabel("F")
        axs[0].set_ylim(0,1)
        axs[0].set_xlim(0,max_time)

        axs[1].plot(times, D.sol(times)[0])
        axs[1].set_ylabel("D")
        axs[1].set_ylim(0,1)
        axs[1].set_xlim(0,max_time)

        axs[2].plot(times, EPSC)
        axs[2].set_ylabel("ESPC")
        axs[2].set_ylim(-1,0)
        axs[2].set_xlim(0,max_time)

        plt.show()
