import numpy as np
from scipy.integrate import solve_ivp
from math import e, floor, factorial

class two_pool(object):
    """class for two pool model
    only works for regular stimuli at r hz"""

    def __init__(self, r, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F):

        EPSC = []
        fastEPSC = []
        slowEPSC = []

        F = 1
        K_F = 1.5
        BG_F = 1e-5
        sat_F = K_F
        facil = BG_F

        fastpool_initial = size_fast
        slowpool_initial = size_slow

        fastpool = fastpool_initial
        slowpool = slowpool_initial

        for i in range(n_pulses):

            facil = facil*e**(-1/(T_F*r))
            F = 1 + sat_F*facil/(K_F + facil)

            fastEPSC.append(fastpool*p_fast*F)
            slowEPSC.append(slowpool*p_slow*F)
            EPSC.append(slowEPSC[-1] + fastEPSC[-1])

            fastpool -= fastpool*p_fast*F
            slowpool -= slowpool*p_slow*F

            fastrecover = (fastpool_initial - fastpool)*(1-e**(-1/(r*T_fast)))
            fastpool += fastrecover
            slowrecover = (slowpool_initial - slowpool)*(1-e**(-1/(r*T_slow)))
            slowpool += slowrecover

            facil += delta_F

        self.norm_val = EPSC[0]
        self.EPSC = np.asarray(EPSC)/EPSC[0]
        self.fastEPSC = np.asarray(fastEPSC)/EPSC[0]
        self.slowEPSC = np.asarray(slowEPSC)/EPSC[0]
        self.final_state = [fastpool, slowpool]

        self.r = r
        self.n_pulses = n_pulses
        self.fastpool_initial = size_fast
        self.slowpool_initial = size_slow
        self.p_fast = p_fast
        self.p_slow = p_slow
        self.T_fast = T_fast
        self.T_slow = T_slow
        self.delta_F = delta_F
        self.T_F = T_F

class maturation(object):

    def __init__(self, r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation):

        n_sites = 1
        delta_t = 1000/r

        #release list preallocation
        n_release = []

        state = np.asarray([0, 0, n_sites, 0], dtype = 'float64') #initial state (empty, immature, mature, facil. sites)
        state = np.reshape(state, (4,1))
        EPSC = []
        EPSC_immature = []
        EPSC_mature = []
        EPSC_facilitated = []

        for i in range(n_pulses):

            n_release.append([state[1]*p_immature, state[2]*p_mature, state[3]*p_facilitated]) #immature, mature, and facilitated release
            EPSC_immature.append(state[1]*p_immature)
            EPSC_mature.append(state[2]*p_mature)
            EPSC_facilitated.append(state[3]*p_facilitated)
            EPSC.append(sum(n_release[-1]))

            state[0] += sum(n_release[-1]) #add released vesicles to empty sites
            state[1:] -= n_release[-1] #account for released vesicles in state change

            n_E_I = state[0]*(1-e**(-1*delta_t/T_refill)) #start in E end in I
            n_E_I_F = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_facilitation)) #start in E end in F after passing through I
            n_E_I_M = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_maturation))
            n_E_I_M_F = state[0]*(1-e**(-1*delta_t/T_refill))*(1-e**(-1*delta_t/T_maturation))*(1-e**(-1*delta_t/T_facilitation))
            n_I_M = state[1]*(1-e**(-1*delta_t/T_maturation))
            n_I_F = state[1]*(1-e**(-1*delta_t/T_facilitation))
            n_I_M_F =  state[1]*(1-e**(-1*delta_t/T_maturation))*(1-e**(-1*delta_t/T_facilitation))
            n_M_F = state[2]*(1-e**(-1*delta_t/T_facilitation))

            state[0] = state[0]*e**(-1*delta_t/T_refill) #after vesicles are released, additional empty sites refill with characteristic time T_refill over delta_t, reducing the number of empty sites
            state[1] += n_E_I - n_I_M - n_E_I_M - n_I_F - n_E_I_F  #the pool begins to refill with immature vesicles, over delta_t, some facilitate, some mature, and a small number do both, some of the old vesicles mature and some facilitate

            state[2] += n_I_M + n_E_I_M - n_M_F - n_I_M_F - n_E_I_M_F #some of the mature pool is replenished by the immature pool that remains after the pulse and some is replenished by the new vesicles refilling pool first, some mature vesicles are lost to facilitation

            state[3] += n_E_I_F + n_E_I_M_F + n_I_F + n_I_M_F + n_M_F #vesicles become facilitated by 5 paths, 2 starting empty, 2 starting immature, and one starting mature

        self.norm_val = EPSC[0]
        self.EPSC = np.asarray(EPSC)/EPSC[0]
        self.EPSC_immature = np.asarray(EPSC_immature)/EPSC[0]
        self.EPSC_mature = np.asarray(EPSC_mature)/EPSC[0]
        self.EPSC_facilitated = np.asarray(EPSC_facilitated)/EPSC[0]

        self.final_state = state

        self.p_immature = p_immature
        self.p_mature = p_mature
        self.p_facilitated = p_facilitated
        self.T_refill = T_refill
        self.T_maturation = T_maturation
        self.T_facilitation = T_facilitation

class maturation_2(object):

    def __init__(self, stimulus_times, step_size, k_rep, k_mat, k_unmat, k_on_1, k_off_1, k_on_2, k_off_2, f, s, k_fuse_basal, Ca_rest, Ca_spike, T_Ca_decay):

        #Ca sensors are modeled as silent for immature vesicles
        #assume that the reserve pool is non depletable
        m = 2 #second sensor cooperativity of 2, god this notation is weird
        n = 5 #syt1 Ca cooperativity of 5
        b = 0.5 #sensor cooperativity

        state = np.zeros((m+1,n+1)) #(i,j) holds the state of the pool with i Ca bound to second sensor and j Ca bound to syt1
        for i in range(m+1):
            for j in range(n+1):
                state[i,j] = (((factorial(n)/factorial(n - j))*(Ca**j)*(k_on_1**j))/(factorial(j)*(b**(j*(j-1)/2))*(k_off_1**j)))*(((factorial(m)/factorial(m - i))*(Ca**i)*(k_on_2**i))/(factorial(i)*(b**(i*(i-1)/2))*(k_off_2**i))) #equation from Kobersmed et al for R(n,m)

        immature =
        Fused = 0
        ###

        def dCa():

            if t in stimulus_times:
                return Ca_spike
            return

        ###

        def partials(t, y):
            immature = y[0]
            
        #Partials with respect to time, number after subscript are number of Ca bound syt 1 and Ca bound second sensor respectively

        def dImmature(): #function for number of immature vesicles

            return(np.sum(k_unmat*state) - k_mat*immature)

        ###

        def dMature_00():

            return(k_rep*(1 - sum(state) - immature) + k_mat*immature - (5*Ca*k_on_1 + 2*Ca*k_on_2 + k_fuse_basal + k_unmat)*state[0,0] + k_off_1*state[1,0] + k_off_2*state[0,1])

        def dMature_10():

            return(-1*(4*Ca*k_on_1 + k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f + k_unmat)*state[1,0] + 5*Ca*k_on_1*state[0,0] + 2*b*k_off_1*state[2,0] + k_off_2*state[1,1])

        def dMature_20():

            return(-1*(3*Ca*k_on_1 + 2*b*k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f**2 + k_unmat)*state[2,0] + 4*Ca*k_on_1*state[1,0] + 3*b*k_off_1*state[3,0] + k_off_2*state[2,1])

        def dMature_30():

            return(-1*(2*Ca*k_on_1 + 3*b**2*k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f**3 + k_unmat)*state[3,0] + 3*Ca*k_on_1*state[2,0] + 4*b**3*k_off_1*state[4,0] + k_off_2*state[3,1])

        def dMature_40():

            return(-1*(Ca*k_on_1 + 4*b**3*k_off_1 + 2*Ca*k_off_2 + k_fuse_basal*f**4 + k_unmat)*state[4,0] + 2*Ca*k_on_1*state[3,0] + 5*b**4*k_off_1*state[5,0] + k_off_2*state[4,1])

        def dMature_50():

            return(-1*(5*b**4*k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f**5 + k_unmat)*state[5,0] + Ca*k_on_1*state[4,0] + k_off_2*state[5,1])

        ###

        def dMature_01():

            return(-1*(5*Ca*k_on_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*s + k_unmat)*state[0,1] + k_off_1*state[1,1] + 2*Ca*k_on_2*state[0,0] + 2*b*k_off_2*state[0,2])

        def dMature_11():

            return(-1*(4*Ca*k_on_1 + k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f*s + k_unmat)*state[1,1] + 5*Ca*k_on_1*state[0,1] + 2*b*k_off_1*state[2,0] + 2*Ca*k_on_2*state[1,0] + 2*b*k_off_2*state[1,2])

        def dMature_21():

            return(-1*(3*Ca*k_on_1 + 2*b*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**2*s + k_unmat)*state[2,1] + 4*Ca*k_on_1*state[1,1] + 3*b**2*k_off_1*state[3,1] + 2*Ca*k_on_2*state[2,0] + 2*b*k_off_2*state[2,2])

        def dMature_31():

            return(-1*(2*Ca*k_on_1 + 3*b**2*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**3*s + k_unmat)*state[3,1] + 3*Ca*k_on_1*state[2,1] + 4*b**3*k_off_1*state[4,1] + 2*Ca*k_on_2*state[3,0] + 2*b*k_off_2*state[3,2])

        def dMature_41():

            return(-1*(Ca*k_on_1 + 4*b**3*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**4*s + k_unmat)*state[4,1] + 2*Ca*k_on_1*state[3,1] + 5*b**3*k_off_1*state[5,1] + 2*Ca*k_on_2*state[4,0] + 2*b*k_off_2*state[4,2])

        def dMature_51():

            return(-1*(5*b**4*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**5*s + k_unmat)*state[5,1] + Ca*k_on_1*state[4,1] + 2*Ca*k_on_2*state[5,0] + 2*b*k_off_2*state[5,2])

        ###

        def dMature_02():

            return(-1*(5*Ca*k_on_1 + 2*b*k_off_2 + k_fuse_basal*s**2 + k_unmat)*state[0,2] + k_off_1*state[1,2] + Ca*k_on_2*state[0,1])

        def dMature_12():

            return(-1*(4*Ca*k_on_1 + k_off_1 + 2*b*k_off_2 + k_fuse_basal*f*s**2 + k_unmat)*state[1,2] + 5*Ca*k_on_1*state[0,2] + 2*b*k_off_1*state[2,2] + Ca*k_on_2*state[1,1])

        def dMature_22():

            return(-1*(3*Ca*k_on_1 + 2*b*k_off_1 + 2*b*k_off_2 + k_fuse_basal*f**2*s**2 + k_unmat)*state[2,2] + 4*Ca*k_on_1*state[1,2] + 3*b**2*k_off_1*state[3,2] + Ca*k_on_2*state[2,1])

        def dMature_32():

            return(-1*(2*Ca*k_on_1 + 3*b**2*k_off_1 + 2*b*k_off_2 + Ca*k_on_2 + k_fuse_basal*f**3*s**2 + k_unmat)*state[3,2] + 3*Ca*k_on_1*state[2,2] + 4*b**3*k_off_1*state[4,2] + Ca*k_on_2*state[3,1])

        def dMature_42():

            return(-1*(Ca*k_on_1 + 4*b**3*k_off_1 + 2*b*k_off_2 + Ca*k_on_2 + k_fuse_basal*f**4*s**2 + k_unmat)*state[4,2] + 2*Ca*k_on_1*state[3,2] + 5*b**3*k_off_1*state[5,2] + Ca*k_on_2*state[4,1])

        def dMature_52():

            return(-1*(5*b**4*k_off_1 + 2*b*k_off_2 + k_fuse_basal*f**5*s**2 + k_unmat)*state[5,2] + Ca*k_on_1*state[4,2] + Ca*k_on_2*state[5,1])

        def dFused():

            return(k_fuse_basal*(state[0,0] + f*state[1,0] + f**2*state[2,0] + f**3*state[3,0] + f**4*state[4,0] + f**5*state[5,0] + s*state[0,1] + f*s*state[1,1] + f**2*s*state[2,1] + f**3*s*state[3,1] + f**4*s*state[4,1] + f**5*s*state[5,1] + s**2*state[0,2] + f*s**2*state[1,2] + f**2*s**2*state[2,2] + f**3*s**2*state[3,2] + f**4*s**2*state[4,2] + f**5*s**2*state[5,2]))





        ts = np.linspace(1,stimulus_times[-1],stimulus_times[-1])

        Ca_1 = Ca_influx*(e**(-1*t/T_Ca_decay) #simnple exponential decay of calcium from max value of Ca_influx
        Ca = np.zeros(stimulus_times[-1]) #functional alpha function

        for i in range(len(stimulus_times)):
            stimulus = stimulus_times[i]
            Ca[stimulus-1:stimulus_times[-1]-1] += alpha_1[0:stimulus_times[-1]-stimulus]) #calculate effect of each stimulus on the Ca level
