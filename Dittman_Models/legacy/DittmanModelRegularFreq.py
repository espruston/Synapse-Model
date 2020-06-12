from math import e

#table 2 parameters
roh = 3.1
F_1 = 0.05
T_F = 100
T_D = 50
k_0 = 2/1000
k_max = 30/1000
K_D = 2

#other parameters
r = 50 #frequency in hz
n_pulses = 20
T_E = 2
delta_F = 1
delta_D = 1

class regular_train(object):
    """Class for stimuli at regular intervals"""

    def __init__(self, n_pulses, r, N_T, roh, F_1, T_F, T_D, K_D, k_0, k_max, K_F, delta_F, delta_D):

        #experiment-dependent
        self.n_pulses = n_pulses #number of pulses to be simulated
        self.r = r #frequency
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
        self.K_F = K_F #affinity of CaX_f for release site
        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #assumed/explicit values
        self.CaX_F = [0] #vector of CaX_F values which is populated by iterate_pulses, CaX_F = 0
        self.F = [F_1] #vector of D values which is populated by iterate_pulses, F_1 determined by user input
        self.CaX_D = [] #vector of CaX_D values which is populated by iterate_pulses
        self.xi = [0] #vector of xi values which is populated by iterate_pulses, xi_1 = 0
        self.D = [1] #vector of D values which is populated by iterate_pulses, D_1 = 1
        self.EPSC_norm = [] #vector of EPSC values normalized to the first pulse which is populated by iterate_pulses

        #steady state values that are used in iteration
        self.CaX_F_ss = self.delta_F*(e**(1/(r*T_F)) - 1)**(-1)
        self.CaX_D_ss = self.delta_D*(1-e**(-1/(r*T_D)))

        #ss values used after iteration
        self.F_ss = F_1 + (1 - F_1)/(1+(self.K_F/self.CaX_F_ss)) #eq 11
        self.xi_ss = ((self.K_D/self.CaX_D_ss + 1)/((self.K_D/self.CaX_D_ss) + e**(-1/(self.r*self.T_D))))**(-1*(self.k_max-self.k_0)*self.T_D) #modified eq 16
        self.D_ss = (1 - e**(-1*self.k_0/self.r)*self.xi_ss)/(1 - (1 - self.F_ss)*e**(-1*self.k_0/self.r)*self.xi_ss) #eq 20
        self.EPSC_norm_ss = self.D_ss*(self.F_ss/F_1) #eq 21

        for i in range(self.n_pulses): #i ranges from 0->npulses-1, so i should be itereated as i+1

            self.CaX_F.append(self.CaX_F_ss*(1-e**(-1*(i)/(self.r*self.T_F))) + self.delta_F) #eq 17

            self.F.append(self.F_1 + (1 - self.F_1)/(1+(self.K_F/self.CaX_F[i+1]))) #eq 10

            self.CaX_D.append(self.CaX_D_ss*(1 - e**(-1*(i)/(self.r*self.T_D))) + self.delta_D) #eq 17

            self.xi.append(((self.K_D/self.CaX_D[i] + 1)/((self.K_D/self.CaX_D[i]) + e**(-1/(self.r*self.T_D))))**(-1*(self.k_max-self.k_0)*self.T_D)) #equation 16

            self.D.append(1 - (1 - (1 - self.F[i]))*(self.D[i])*e**(self.k_0/self.r)*self.xi[i+1]) #equation 15

            self.EPSC_norm.append(self.D[i+1]*(self.F[i+1]/self.F_1)) #eq 19


        self.CaX_D.append(self.CaX_D_ss*(1-e**(-1*(self.n_pulses)/(self.r*self.T_D)))) #calculate value after last pulse for the last pulse

        def plot_EPSCs(self, ts):
            """ts should be an evenly spaced, strictly increasing vector, starting at 0, containing the times which are to be plotted on the x axis"""



            alpha = (t*e/T_E)*e**(-t/T_E)

            return

class poisson_train(object):
    """Class for stimuli according to a poisson train"""

    def __init__(self, delta_ts, roh, F_1, T_F, T_D, K_D, k_0, k_max, delta_ts, N_T, K_F, delta_F, delta_D):

        #experiment-dependent
        self.delta_ts = delta_ts #time steps between pulses, first pulse occurs after delta_ts[0] msec
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
        self.K_F = K_F #affinity of CaX_f for release site
        self.delta_F = delta_F #amount of CaX_F increase as a result of stimulus
        self.delta_D = delta_D #amount of CaX_D increase as a result of stimulus

        #assumed/explicit values
        self.CaX_F = [0] #vector of CaX_F values which is populated by iterate_pulses, CaX_F = 0
        self.F = [F_1] #vector of D values which is populated by iterate_pulses, F_1 determined by user input
        self.CaX_D = [] #vector of CaX_D values which is populated by iterate_pulses
        self.xi = [0] #vector of xi values which is populated by iterate_pulses, xi_1 = 0
        self.D = [1] #vector of D values which is populated by iterate_pulses, D_1 = 1
        self.EPSC = [] #vector of EPSC values normalized to the first pulse which is populated by iterate_pulses

        for i in range(len(self.delta_ts)):

            self.CaX_F.append(self.CaX_F[i]*e**(-1*self.delta_ts[i+1]/self.T_F) + self.delta_F)

            self.F.append(self.F_1 + (1-self.F_1)/(1+(self.K_F/self.CaX_F[i+1])))

            self.CaX_D.append(self.CaX_D[i]*e**(-1*self.delta_ts[i+1]/self.T_D) + self.delta_D)

            self.xi.append(((self.K_D/self.CaX_D[i] + 1)/((self.K_D/self.CaX_D[i]) + e**(-1/(self.r*self.T_D))))**(-1*(self.k_max-self.k_0)*self.T_D)) #equation 16

            self.D.append(1 - (1 - (1 - self.F[i]))*(self.D[i])*e**(self.k_0/self.r)*self.xi[i+1])

            self.EPSC.append(self.N_T*self.D[i+1]*self.F[i+1])

        self.CaX_D.append(self.CaX_D[len(self.delta_ts)]*e**(-1*self.delta_ts[len(self.delta_ts)+1]/self.T_D) + self.delta_D)


    def plot_EPSCs(self, ts):
        """ts should be an evenly spaced, strictly increasing vector, starting at 0, containing the times which are to be plotted on the x axis"""



        alpha = (t*e/T_E)*e**(-t/T_E)

        return
