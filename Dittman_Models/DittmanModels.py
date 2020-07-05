import numpy as np
from scipy.integrate import solve_ivp
from math import e

class regular_train(object):
    """Class for stimuli at regular intervals"""

    def __init__(self, stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):

        #experiment-dependent
        self.stimulus_times = stimulus_times
        self.n_pulses = len(stimulus_times) #number of pulses to be simulated
        self.r = 1/np.diff(stimulus_times)[0] #Hz in msec
        self.max_time = max_time
        self.N_T = N_T #number of total release sites

        #table 2 values
        self.roh = roh
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0/1000
        self.k_max = k_max/1000
        self.K_D = K_D

        #testing parameters
        self.K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*roh - F_1) - 1) #eq 7, affinity of CaX_f for release site
        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #assumed/explicit values,
        self.CaX_F = [] #vector of CaX_F values, CaX_F[0] = 0
        self.F = [F_1] #vector of F values
        self.CaX_D = [0] #vector of CaX_D values, CaX_D[0] is basically calculated twice (once here, and once in the loop) due to eq 17
        self.xi = [] #vector of xi values, xi_1 = 1
        self.D = [1] #vector of D values, D_1 = 1
        self.EPSC = [] #vector of EPSC values

        #steady state values that are used in iteration
        if 1/(self.r*self.T_F) > 709: #prevent overflow
            self.CaX_F_ss = 0
        else:
            self.CaX_F_ss = self.delta_F*((e**(1/(self.r*self.T_F)) - 1)**(-1))
        if 1/(self.r*self.T_D) > 709:
            self.CaX_D_ss = 0
        else:
            self.CaX_D_ss = self.delta_D*((1 - e**(-1/(self.r*self.T_D)))**(-1))

        #ss values used after iteration
        if self.CaX_F_ss == 0:
            self.F_ss = F_1
        else:
            self.F_ss = F_1 + (1 - F_1)/(1+(self.K_F/self.CaX_F_ss)) #eq 11

        if self.CaX_D_ss == 0:
            self.xi_ss = 1
        else:
            self.xi_ss = ((self.K_D/self.CaX_D_ss + 1)/((self.K_D/self.CaX_D_ss) + e**(-1/(self.r*self.T_D))))**(-1*(self.k_max-self.k_0)*self.T_D) #modified eq 16
        self.D_ss = (1 - e**(-1*self.k_0/self.r)*self.xi_ss)/(1 - (1 - self.F_ss)*e**(-1*self.k_0/self.r)*self.xi_ss) #eq 20
        self.EPSC_norm_ss = self.D_ss*(self.F_ss/F_1) #eq 21

        for i in range(self.n_pulses): #i ranges from 0->npulses-1

            self.CaX_F.append(self.CaX_F_ss*(1-e**(-1*(i)/(self.r*self.T_F)))) #eq 8 amount of CaX_F just before current stimulus

            if self.CaX_F[-1] == 0:
                self.F.append(self.F_1)
            else:
                self.F.append(self.F_1 + (1 - self.F_1)/(1+(self.K_F/self.CaX_F[-1]))) #eq 10

            self.CaX_D.append(self.CaX_D_ss*(1 - e**(-1*(i)/(self.r*self.T_D)))) #eq 17, actually CaX_D_(i-1) in the Dittman paper, amount of CaX_D just after previous pulse

            #in calculating xi at pulse n = i, CaX_D[-1] defines CaX_D after pulse i (previous pulse)
            if self.CaX_D[-1] == 0:
                self.xi.append(1)
            else:
                self.xi.append(((self.K_D/self.CaX_D[-1] + 1)/((self.K_D/self.CaX_D[-1]) + e**(-1/(self.r*self.T_D))))**(-1*(self.k_max-self.k_0)*self.T_D)) #equation 16

            self.D.append(1 - (1 - (1 - self.F[-2])*(self.D[-1]))*e**(-1*self.k_0/self.r)*self.xi[-1]) #equation 15

            self.EPSC.append(self.N_T*self.D[-1]*self.F[-1]) #eq 19, first pulse at stimulus_times[0]


        self.CaX_D.append(self.CaX_D_ss*(1 - e**(-1*(self.n_pulses-1)/(self.r*self.T_D)))) #calculate value after last pulse for the last pulse

        self.times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        alpha_1 = (self.times*e/T_E)*e**(-1*self.times/T_E) #reference alpha function
        self.alpha = np.zeros(self.max_time) #functional alpha function

        for i in range(len(self.stimulus_times)):
            stimulus = self.stimulus_times[i]
            self.alpha[stimulus-1:self.max_time-1] = self.alpha[stimulus-1:self.max_time-1] + (self.F[i+1] * self.D[i+1] * alpha_1[0:self.max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func = self.alpha * self.N_T


class poisson_train(object):
    """Class for stimuli according to a poisson train"""

    def __init__(self, stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):

        #experiment-dependent
        self.stimulus_times = stimulus_times #time steps between pulses, first pulse occurs after delta_ts[0] msec
        self.delta_ts = np.concatenate(([stimulus_times[0]], np.diff(stimulus_times)))
        self.N_T = N_T #number of total release sites
        self.max_time = max_time

        #table 2 values
        self.roh = roh
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0/1000
        self.k_max = k_max/1000
        self.K_D = K_D

        #testing parameters
        self.K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*roh - F_1) - 1) #affinity of CaX_f for release site
        if self.K_F == 0:
            self.K_F = 1e30

        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #assumed/explicit values
        self.CaX_F = [0] #vector of CaX_F values, CaX_F = 0
        self.F = [F_1] #vector of D values, F_1 determined by user input
        self.CaX_D = [0] #vector of CaX_D values
        self.xi = [1] #vector of xi values, xi_1 = 1
        self.D = [1] #vector of D values, D_1 = 1
        self.EPSC = [] #vector of EPSC values normalized to the first pulse

        for i in range(len(self.stimulus_times)):

            self.CaX_F.append(self.CaX_F[-1]*e**(-1*self.delta_ts[i]/self.T_F) + self.delta_F)

            if self.CaX_F[-2] == 0:
                self.F.append(self.F_1)
            else:
                self.F.append(self.F_1 + (1-self.F_1)/(1+(self.K_F/self.CaX_F[-2])))

            self.CaX_D.append(self.CaX_D[-1]*e**(-1*self.delta_ts[i]/self.T_D) + self.delta_D)

            if self.CaX_D[-2] == 0:
                self.xi.append(1)
            else:
                self.xi.append(((self.K_D/self.CaX_D[-2] + 1)/((self.K_D/self.CaX_D[-2]) + e**(-1/((1/self.delta_ts[i])*self.T_D))))**(-1*(self.k_max - self.k_0)*self.T_D)) #equation 16, 1/delta_t is the frequency in msec for defining the interval between pulse i and i-1

            self.D.append(1 - (1 - (1 - self.F[-2])*(self.D[-1]))*e**(-1*self.k_0/(1/self.delta_ts[i]))*self.xi[-1])

            self.EPSC.append(self.N_T*self.D[-1]*self.F[-1])


        self.CaX_D.append(self.CaX_D[-1]*e**(-1*self.delta_ts[-1]/self.T_D) + self.delta_D)

        self.times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        alpha_1 = (self.times*e/T_E)*e**(-1*self.times/T_E) #reference alpha function
        self.alpha = np.zeros(self.max_time) #functional alpha function

        for i in range(len(self.stimulus_times)):
            stimulus = self.stimulus_times[i]
            self.alpha[stimulus-1:self.max_time-1] = self.alpha[stimulus-1:self.max_time-1] + (self.F[i+1] * self.D[i+1] * alpha_1[0:self.max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func = self.alpha * self.N_T


class DittmanRK1(object):

    def __init__(self, stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):
        #experiment-dependent
        self.stimulus_times = stimulus_times #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector!!
        self.max_time = max_time
        self.N_T = N_T #number of total release sites

        #table 2 values
        self.roh = roh
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0/1000
        self.k_max = k_max/1000
        self.K_D = K_D

        #testing parameters
        self.K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*roh - F_1) - 1) #affinity of CaX_f for release site
        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #set or calculated values
        self.times = np.linspace(1,self.max_time,self.max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        self.stimuli = np.zeros(max_time)
        for stimulus in self.stimulus_times:
            self.stimuli[stimulus] = 1 #stiumli vector now is 1 for times when a stimulus occurs and 0 if it doesn't

        #populate CaX_F
        self.CaX_F = np.zeros(self.max_time+1) #preallocate vector of CaX_F values from t=0->max_time
        #self.CaX_F = self.CaX_F + 1e-30 #avoid /0 errors
        for t in self.times:
            dCaX_F_dt = -1*self.CaX_F[t-1]/self.T_F + self.delta_F*self.stimuli[t-1] #exponential decay with characteristic time T_F after a jump of delta_F after an action potential (Eq. 3)
            self.CaX_F[t] = self.CaX_F[t-1] + dCaX_F_dt #CaX_F[0] should be 0 so that no facilitation occurs before an action potential

        self.F = np.zeros(self.max_time+1) #preallocate vector of F values from t=0->max_time
        self.F[0] = self.F_1 #F = F_1 until a stimulus occurs
        for t in self.times:
            if self.CaX_F[t] == 0:
                self.F[t] = self.F_1
            else:
                self.F[t] = self.F_1 + (1 - self.F_1)/(1 + self.K_F/self.CaX_F[t]) #(Eq. 2a)

        self.CaX_D = np.zeros(self.max_time+1) #preallocate vector of CaX_D values from t=0->max_time
        for t in self.times:
            dCaX_D_dt = -1*self.CaX_D[t-1]/self.T_D + self.delta_D*self.stimuli[t-1] #exponential decay with characteristic time T_D after a jump of delta_D after an action potential (Eq. 12)
            self.CaX_D[t] = self.CaX_D[t-1] + dCaX_D_dt #CaX_D[0] should be 0 so that np Ca dependent recovery occurs before an action potential

        self.k_recov = np.zeros(self.max_time+1) #preallocate values of k_recov
        self.k_recov[0] = self.k_0
        for t in self.times:
            if self.CaX_D[t] ==0:
                self.k_recov[t] = self.k_0
            else:
                self.k_recov[t] = (self.k_max - self.k_0)/(1 + self.K_D/self.CaX_D[t]) + self.k_0

        self.D = np.ones(self.max_time+1) #preallocate vector of D values from t=0->max_time (D=1 @ t=0)
        for t in self.times:
            dD_dt = (1 - self.D[t-1])*self.k_recov[t-1] - self.D[t-1]*self.F[t-1]*self.stimuli[t-1] #scale k_recov to msec
            #dD_dt = (1 - D[t-1])*k_recov/1000 - D[t-1]*F[t-1]*stimuli[t-1]
            self.D[t] = self.D[t-1] + dD_dt

        self.EPSC = self.F[(self.stimulus_times)-1]*self.D[(self.stimulus_times)-1]

        alpha_1 = (self.times*e/T_E)*e**(-1*self.times/T_E) #reference alpha function
        self.alpha = np.zeros(self.max_time) #functional alpha function

        for stimulus in self.stimulus_times:
            self.alpha[stimulus-1:self.max_time-1] = self.alpha[stimulus-1:self.max_time-1] + (self.F[stimulus-1] * self.D[stimulus-1] * alpha_1[0:self.max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func =  self.alpha * self.N_T


class DittmanRK45(object):

    def __init__(self, stimulus_times, max_time, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_F, delta_D, T_E):
        #experiment-dependent
        self.stimulus_times = stimulus_times #vector containing times at which stimuli occur (in msec) !!make sure these values exist in times vector!!
        self.max_time = max_time
        self.N_T = N_T #number of total release sites

        #table 2 values
        self.roh = roh
        self.F_1 = F_1
        self.T_F = T_F
        self.T_D = T_D
        self.k_0 = k_0/1000
        self.k_max = k_max/1000
        self.K_D = K_D

        #testing parameters
        self.K_F = delta_F*((1 - F_1)/((F_1/(1 - F_1))*roh - F_1) - 1) #affinity of CaX_f for release site
        if self.K_F == 0:
            self.K_F = 1e10

        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #set or calculated values
        self.max_step_size = 1
        self.times = np.linspace(1,max_time,max_time,dtype=np.int32) #vector of length max_time denoting times from 1->max_time msec with step size 1

        #partials to be solved
        def dCaX_F_dt(t, CaX_F):

            return -1*CaX_F/self.T_F + self.delta_F*(int(t) in self.stimulus_times) #exponential decay with characteristic time T_F after a jump of delta_F after an action potential (Eq. 3)

        def dCaX_D_dt(t, CaX_D):

            return -1*CaX_D/self.T_D + self.delta_D*(int(t) in self.stimulus_times) #exponential decay with characteristic time T_D after a jump of delta_D after an action potential (Eq. 12)

        def dD_dt(t, D):

            if self.CaX_D.sol(t) == 0:
                k_recov_temp = self.k_0
            else:
                k_recov_temp = (self.k_max - self.k_0)/(1 + self.K_D/self.CaX_D.sol(t)[0]) + self.k_0 #eq 14

            return (1 - D)*k_recov_temp - D*self.F[int(t)-1]*(int(t) in self.stimulus_times) #eq 13

        #values solved for in ivp
        self.CaX_F = solve_ivp(dCaX_F_dt, [0, self.max_time], [0], t_eval = self.stimulus_times, dense_output = True, first_step = self.stimulus_times[0], max_step = self.max_step_size) #solve dCaX_F_dt using RK45 starting from first stimulus w max step size 1
        self.CaX_D = solve_ivp(dCaX_D_dt, [0, max_time], [0], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = self.max_step_size) #solve dCaX_D_dt using RK45 starting from first stimulus w max step size 1

        self.F = self.F_1 + (1 - self.F_1)/(1 + self.K_F/(self.CaX_F.sol(self.times)[0] + 1e-30*(self.CaX_F.sol(self.times)[0] == 0))) #equation 2

        self.k_recov = (self.k_max - self.k_0)/(1 + self.K_D/(self.CaX_D.sol(self.times)[0] + 1e-30*(self.CaX_D.sol(self.times)[0] == 0))) + self.k_0
        self.D_object = solve_ivp(dD_dt, [0, max_time], [1], t_eval = stimulus_times, dense_output = True, first_step = stimulus_times[0], max_step = self.max_step_size) #solve dD_dt using RK45 starting from first stimulus w max step size 1
        self.D = self.D_object.sol(self.times)[0]

        self.EPSC = self.F[(self.stimulus_times)-1]*self.D[(self.stimulus_times)-1]

        alpha_1 = (self.times*e/T_E)*e**(-1*self.times/T_E) #reference alpha function
        self.alpha = np.zeros(self.max_time) #functional alpha function

        for stimulus in self.stimulus_times:
            self.alpha[stimulus-1:self.max_time-1] = self.alpha[stimulus-1:self.max_time-1] + (self.F[stimulus-1] * self.D[stimulus-1] * alpha_1[0:self.max_time-stimulus]) #calculate effect of each stimulus on the alpha function and sum them

        self.EPSC_func =  self.alpha * self.N_T
