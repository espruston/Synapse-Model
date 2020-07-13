import numpy as np
from scipy.integrate import solve_ivp
from math import e

class regular_train(object):
    """Class for stimuli at regular intervals"""

    def __init__(self, stimulus_times, max_time, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):

        #experiment-dependent
        n_pulses = len(stimulus_times) #number of pulses to be simulated
        r = 1/np.diff(stimulus_times)[0] #Hz in msec
        #table 2 values

        k_0 = k_0/1000
        k_max = k_max/1000

        #testing parameters
        K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*rho - F_1) - 1) #eq 7, affinity of CaX_f for release site

        #assumed/explicit values,
        CaX_F = [] #vector of CaX_F values, CaX_F[0] = 0
        F = [F_1] #vector of F values
        CaX_D = [0] #vector of CaX_D values, CaX_D[0] is basically calculated twice (once here, and once in the loop) due to eq 17
        xi = [] #vector of xi values, xi_1 = 1
        D = [1] #vector of D values, D_1 = 1
        EPSC = [] #vector of EPSC values

        #steady state values that are used in iteration
        if 1/(r*T_F) > 709: #prevent overflow
            CaX_F_ss = 0
        else:
            CaX_F_ss = delta_F*((e**(1/(r*T_F)) - 1)**(-1))
        if 1/(r*T_D) > 709:
            CaX_D_ss = 0
        else:
            CaX_D_ss = delta_D*((1 - e**(-1/(r*T_D)))**(-1))

        #ss values used after iteration
        if CaX_F_ss == 0:
            F_ss = F_1
        else:
            F_ss = F_1 + (1 - F_1)/(1+(K_F/CaX_F_ss)) #eq 11

        if CaX_D_ss == 0:
            xi_ss = 1
        else:
            try:
                xi_ss = ((K_D/CaX_D_ss + 1)/((K_D/CaX_D_ss) + e**(-1/(r*T_D))))**(-1*(k_max-k_0)*T_D) #modified eq 16
            except:
                print(K_D, CaX_D_ss, r, T_D, k_max, k_0)
        D_ss = (1 - e**(-1*k_0/r)*xi_ss)/(1 - (1 - F_ss)*e**(-1*k_0/r)*xi_ss) #eq 20
        EPSC_norm_ss = D_ss*(F_ss/F_1) #eq 21

        for i in range(n_pulses): #i ranges from 0->npulses-1
            EPSC.append(N_T*D[-1]*F[-1]) #eq 19, first pulse at stimulus_times[0]

            CaX_F.append(CaX_F_ss*(1-e**(-1*(i)/(r*T_F)))) #eq 8 amount of CaX_F just before current stimulus

            if CaX_F[-1] == 0:
                F.append(F_1)
            else:
                F.append(F_1 + (1 - F_1)/(1+(K_F/CaX_F[-1]))) #eq 10

            CaX_D.append(CaX_D_ss*(1 - e**(-1*(i)/(r*T_D)))) #eq 17, actually CaX_D_(i-1) in the Dittman paper, amount of CaX_D just after previous pulse

            #in calculating xi at pulse n = i, CaX_D[-1] defines CaX_D after pulse i (previous pulse)
            if CaX_D[-1] == 0:
                xi.append(1)
            else:
                xi.append(((K_D/CaX_D[-1] + 1)/((K_D/CaX_D[-1]) + e**(-1/(r*T_D))))**(-1*(k_max-k_0)*T_D)) #equation 16

            D.append(1 - (1 - (1 - F[-2])*(D[-1]))*e**(-1*k_0/r)*xi[-1]) #equation 15

        CaX_D.append(CaX_D_ss*(1 - e**(-1*(n_pulses-1)/(r*T_D)))) #calculate value after last pulse for the last pulse

        times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        # alpha_1 = (times*e/T_E)*e**(-1*times/T_E) #reference alpha function
        # alpha = np.zeros(max_time) #functional alpha function
        #
        # for i in range(len(stimulus_times)):
        #     stimulus = stimulus_times[i]
        #     alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (F[i+1] * D[i+1] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them
        #
        # self.EPSC_func = alpha * N_T
        self.times = times

        self.stimulus_times = stimulus_times
        self.n_pulses = n_pulses #number of pulses to be simulated
        self.r = r #Hz in msec
        self.max_time = max_time
        self.N_T = N_T #number of total release sites

        #table 2 values
        self.rho = rho
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0
        self.k_max = k_max
        self.K_D = K_D

        #testing parameters
        self.K_F = K_F #eq 7, affinity of CaX_f for release site
        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #assumed/explicit values,
        self.CaX_F_ss = CaX_F_ss
        self.CaX_F = CaX_F #vector of CaX_F values, CaX_F[0] = 0
        self.D = D
        self.F = F #vector of F values
        self.CaX_D_ss = CaX_D_ss
        self.CaX_D = CaX_D #vector of CaX_D values, CaX_D[0] is basically calculated twice (once here, and once in the loop) due to eq 17
        self.xi_ss = xi_ss
        self.xi = xi #vector of xi values, xi_1 = 1
        self.EPSC = np.asarray(EPSC)/EPSC[0]

class poisson_train(object):
    """Class for stimuli according to a poisson train"""

    def __init__(self, stimulus_times, max_time, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):

        #experiment-dependent
        delta_ts = np.concatenate(([stimulus_times[0]], np.diff(stimulus_times)))

        k_0 = k_0/1000
        k_max = k_max/1000

        #testing parameters
        K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*rho - F_1) - 1) #affinity of CaX_f for release site
        if K_F == 0:
            K_F = 1e30

        #assumed/explicit values
        CaX_F = [0] #vector of CaX_F values, CaX_F = 0
        F = [F_1] #vector of D values, F_1 determined by user input
        CaX_D = [0] #vector of CaX_D values
        xi = [1] #vector of xi values, xi_1 = 1
        D = [1] #vector of D values, D_1 = 1
        EPSC = [] #vector of EPSC values normalized to the first pulse

        for i in range(len(stimulus_times)):
            EPSC.append(N_T*D[-1]*F[-1])

            CaX_F.append(CaX_F[-1]*e**(-1*delta_ts[i]/T_F) + delta_F)

            if CaX_F[-2] == 0:
                F.append(F_1)
            else:
                F.append(F_1 + (1-F_1)/(1+(K_F/CaX_F[-2])))

            CaX_D.append(CaX_D[-1]*e**(-1*delta_ts[i]/T_D) + delta_D)

            if CaX_D[-2] == 0:
                xi.append(1)
            else:
                xi.append(((K_D/CaX_D[-2] + 1)/((K_D/CaX_D[-2]) + e**(-1/((1/delta_ts[i])*T_D))))**(-1*(k_max - k_0)*T_D)) #equation 16, 1/delta_t is the frequency in msec for defining the interval between pulse i and i-1

            D.append(1 - (1 - (1 - F[-2])*(D[-1]))*e**(-1*k_0/(1/delta_ts[i]))*xi[-1])


        CaX_D.append(CaX_D[-1]*e**(-1*delta_ts[-1]/T_D) + delta_D)

        times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        alpha_1 = (times*e/T_E)*e**(-1*times/T_E) #reference alpha function
        alpha = np.zeros(max_time) #functional alpha function

        for i in range(len(stimulus_times)):
            stimulus = stimulus_times[i]
            alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (F[i+1] * D[i+1] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func = alpha * N_T
        self.stimulus_times = stimulus_times #time steps between pulses, first pulse occurs after delta_ts[0] msec
        self.delta_ts = delta_ts
        self.N_T = N_T #number of total release sites
        self.max_time = max_time

        #table 2 values
        self.rho = rho
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0
        self.k_max = k_max
        self.K_D = K_D

        #testing parameters
        self.K_F = K_F #affinity of CaX_f for release site

        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #assumed/explicit values
        self.CaX_F = CaX_F #vector of CaX_F values, CaX_F = 0
        self.F = F #vector of D values, F_1 determined by user input
        self.CaX_D = CaX_D #vector of CaX_D values
        self.xi = xi #vector of xi values, xi_1 = 1
        self.D = D #vector of D values, D_1 = 1
        self.EPSC = np.asarray(EPSC)/EPSC[0] #vector of EPSC values normalized to the first pulse


class DittmanRK1(object):

    def __init__(self, stimulus_times, max_time, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):

        k_0 = k_0/1000
        k_max = k_max/1000

        #testing parameters
        K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*rho - F_1) - 1) #affinity of CaX_f for release site

        #set or calculated values
        times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        stimuli = np.zeros(max_time)
        for stimulus in stimulus_times:
            stimuli[stimulus] = 1 #stiumli vector now is 1 for times when a stimulus occurs and 0 if it doesn't

        #populate CaX_F
        CaX_F = np.zeros(max_time+1) #preallocate vector of CaX_F values from t=0->max_time
        #CaX_F = CaX_F + 1e-30 #avoid /0 errors
        for t in times:
            dCaX_F_dt = -1*CaX_F[t-1]/T_F + delta_F*stimuli[t-1] #exponential decay with characteristic time T_F after a jump of delta_F after an action potential (Eq. 3)
            CaX_F[t] = CaX_F[t-1] + dCaX_F_dt #CaX_F[0] should be 0 so that no facilitation occurs before an action potential

        F = np.zeros(max_time+1) #preallocate vector of F values from t=0->max_time
        F[0] = F_1 #F = F_1 until a stimulus occurs
        for t in times:
            if CaX_F[t] == 0:
                F[t] = F_1
            else:
                F[t] = F_1 + (1 - F_1)/(1 + K_F/CaX_F[t]) #(Eq. 2a)

        CaX_D = np.zeros(max_time+1) #preallocate vector of CaX_D values from t=0->max_time
        for t in times:
            dCaX_D_dt = -1*CaX_D[t-1]/T_D + delta_D*stimuli[t-1] #exponential decay with characteristic time T_D after a jump of delta_D after an action potential (Eq. 12)
            CaX_D[t] = CaX_D[t-1] + dCaX_D_dt #CaX_D[0] should be 0 so that np Ca dependent recovery occurs before an action potential

        k_recov = np.zeros(max_time+1) #preallocate values of k_recov
        k_recov[0] = k_0
        for t in times:
            if CaX_D[t] ==0:
                k_recov[t] = k_0
            else:
                k_recov[t] = (k_max - k_0)/(1 + K_D/CaX_D[t]) + k_0

        D = np.ones(max_time+1) #preallocate vector of D values from t=0->max_time (D=1 @ t=0)
        for t in times:
            dD_dt = (1 - D[t-1])*k_recov[t-1] - D[t-1]*F[t-1]*stimuli[t-1] #scale k_recov to msec
            #dD_dt = (1 - D[t-1])*k_recov/1000 - D[t-1]*F[t-1]*stimuli[t-1]
            D[t] = D[t-1] + dD_dt

        EPSC = F[stimulus_times-1]*D[stimulus_times-1]

        alpha_1 = (times*e/T_E)*e**(-1*times/T_E) #reference alpha function
        alpha = np.zeros(max_time) #functional alpha function

        for stimulus in stimulus_times:
            alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (F[stimulus-1] * D[stimulus-1] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func =  alpha * N_T
        self.EPSC = EPSC/EPSC[0]
        self.stimulus_times = stimulus_times #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector!!
        self.max_time = max_time
        self.N_T = N_T #number of total release sites

        #table 2 values
        self.rho = rho
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0
        self.k_max = k_max
        self.k_recov = k_recov
        self.K_D = K_D

        self.F = F
        self.D = D
        self.CaX_F = CaX_F
        self.CaX_D = CaX_D
        #testing parameters
        self.K_F = K_F #affinity of CaX_f for release site
        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        self.stimuli = stimuli
        #set or calculated values
        self.times = times #vector of length max_time denoting times from 1->max_time msec with step size 1


class DittmanRK45(object):

    def __init__(self, stimulus_times, max_time, N_T, rho, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):

        k_0 = k_0/1000
        k_max = k_max/1000

        K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*rho - F_1) - 1) #affinity of CaX_f for release site
        if K_F == 0:
            K_F = 1e10

        max_step_size = 1
        times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        #partials to be solved
        def dCaX_F_dt(t, CaX_F):

            return -1*CaX_F/T_F + delta_F*(int(t) in stimulus_times) #exponential decay with characteristic time T_F after a jump of delta_F after an action potential (Eq. 3)

        def dCaX_D_dt(t, CaX_D):

            return -1*CaX_D/T_D + delta_D*(int(t) in stimulus_times) #exponential decay with characteristic time T_D after a jump of delta_D after an action potential (Eq. 12)

        def dD_dt(t, D):

            if CaX_D.sol(t) == 0:
                k_recov_temp = k_0
            else:
                k_recov_temp = (k_max - k_0)/(1 + K_D/CaX_D.sol(t)[0]) + k_0 #eq 14

            return (1 - D)*k_recov_temp - D*F[int(t)-1]*(int(t) in stimulus_times) #eq 13

        #values solved for in ivp
        CaX_F = solve_ivp(dCaX_F_dt, [0, max_time], [0], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = max_step_size) #solve dCaX_F_dt using RK45 starting from first stimulus w max step size 1
        CaX_D = solve_ivp(dCaX_D_dt, [0, max_time], [0], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = max_step_size) #solve dCaX_D_dt using RK45 starting from first stimulus w max step size 1

        F = F_1 + (1 - F_1)/(1 + K_F/(CaX_F.sol(times)[0] + 1e-30*(CaX_F.sol(times)[0] == 0))) #equation 2

        k_recov = (k_max - k_0)/(1 + K_D/(CaX_D.sol(times)[0] + 1e-30*(CaX_D.sol(times)[0] == 0))) + k_0
        D_object = solve_ivp(dD_dt, [0, max_time], [1], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = max_step_size) #solve dD_dt using RK45 starting from first stimulus w max step size 1
        D = D_object.sol(times)[0]

        EPSC = F[(stimulus_times)-1]*D[(stimulus_times)-1]

        alpha_1 = (times*e/T_E)*e**(-1*times/T_E) #reference alpha function
        alpha = np.zeros(max_time) #functional alpha function

        for stimulus in stimulus_times:
            alpha[stimulus-1:max_time-1] = alpha[stimulus-1:max_time-1] + (F[stimulus-1] * D[stimulus-1] * alpha_1[0:max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func =  alpha * N_T
        self.CaX_F = CaX_F
        self.CaX_D = CaX_D
        self.F = F
        self.D = D
        self.EPSC = np.asarray(EPSC)/EPSC[0]

        self.stimulus_times = stimulus_times #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector!!
        self.max_time = max_time
        self.N_T = N_T #number of total release sites

        #table 2 values
        self.rho = rho
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0
        self.k_max = k_max
        self.k_recov = k_recov
        self.K_D = K_D

        #testing parameters
        self.K_F = K_F #affinity of CaX_f for release site

        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #set or calculated values
        self.max_step_size = max_step_size
        self.times = times #vector of length max_time denoting times from 1->max_time msec with step size 1
