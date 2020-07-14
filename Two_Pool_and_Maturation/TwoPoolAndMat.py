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

    def __init__(self, stimulus_times, k_rep, k_mat, k_unmat, k_on_1, k_off_1, k_on_2, k_off_2, f, s, k_fuse_basal, Ca_rest, Ca_spike, T_Ca_decay):

        #Ca sensors are modeled as silent for immature vesicles
        #assume that the reserve pool is non depletable
        m = 2 #second sensor cooperativity of 2, god this notation is weird
        n = 5 #syt1 Ca cooperativity of 5
        b = 0.5 #sensor cooperativity
        n_sites = 1

        state_ini = np.zeros((m+1)*(n+1)+3) #vectorized for solve_ivp usability
        c = 0 #counter for vector index
        for i in range(m+1):
            for j in range(n+1):
                state_ini[c] = (((factorial(n)/factorial(n - j))*(Ca_rest**j)*(k_on_1**j))/(factorial(j)*(b**(j*(j-1)/2))*(k_off_1**j)))*(((factorial(m)/factorial(m - i))*(Ca_rest**i)*(k_on_2**i))/(factorial(i)*(b**(i*(i-1)/2))*(k_off_2**i))) #equation from Kobersmed et al for R(n,m)
                c += 1

        state_ini[-3] = n_sites - sum(state_ini[0:-3]) #number of immature vesicles
        state_ini[-2] = 0 #number of fused vesicles
        state_ini[-1] = Ca_rest
        ###

        def dCa(t, Ca):

            return (-1/T_Ca_decay)*Ca + Ca_spike*(t in stimulus_times) #spike calcium by Ca_spike at each pulse, otherwise, decay to 0 with characteristic time T_Ca_decay

        ###

        #Partials with respect to time, number after subscript are number of Ca bound syt 1 and Ca bound second sensor respectively

        def dImmature(state): #function for number of immature vesicles

            return(k_unmat*np.sum(state[0:-2]) - k_mat*state[-2])

        ###

        def dMature_00(state, Ca):

            return(k_rep*(n_sites - sum(state[0:-1])) + k_mat*state[-2] - (5*Ca*k_on_1 + 2*Ca*k_on_2 + k_fuse_basal + k_unmat)*state[0] + k_off_1*state[1] + k_off_2*state[6])

        def dMature_10(state, Ca):

            return(-1*(4*Ca*k_on_1 + k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f + k_unmat)*state[1] + 5*Ca*k_on_1*state[0] + 2*b*k_off_1*state[2] + k_off_2*state[7])

        def dMature_20(state, Ca):

            return(-1*(3*Ca*k_on_1 + 2*b*k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f**2 + k_unmat)*state[2] + 4*Ca*k_on_1*state[1] + 3*b*k_off_1*state[3] + k_off_2*state[8])

        def dMature_30(state, Ca):

            return(-1*(2*Ca*k_on_1 + 3*b**2*k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f**3 + k_unmat)*state[3] + 3*Ca*k_on_1*state[2] + 4*b**3*k_off_1*state[4] + k_off_2*state[9])

        def dMature_40(state, Ca):

            return(-1*(Ca*k_on_1 + 4*b**3*k_off_1 + 2*Ca*k_off_2 + k_fuse_basal*f**4 + k_unmat)*state[4] + 2*Ca*k_on_1*state[3] + 5*b**4*k_off_1*state[5] + k_off_2*state[10])

        def dMature_50(state, Ca):

            return(-1*(5*b**4*k_off_1 + 2*Ca*k_on_2 + k_fuse_basal*f**5 + k_unmat)*state[5] + Ca*k_on_1*state[4] + k_off_2*state[11])

        ###

        def dMature_01(state, Ca):

            return(-1*(5*Ca*k_on_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*s + k_unmat)*state[6] + k_off_1*state[7] + 2*Ca*k_on_2*state[0] + 2*b*k_off_2*state[12])

        def dMature_11(state, Ca):

            return(-1*(4*Ca*k_on_1 + k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f*s + k_unmat)*state[7] + 5*Ca*k_on_1*state[6] + 2*b*k_off_1*state[2] + 2*Ca*k_on_2*state[1] + 2*b*k_off_2*state[13])

        def dMature_21(state, Ca):

            return(-1*(3*Ca*k_on_1 + 2*b*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**2*s + k_unmat)*state[8] + 4*Ca*k_on_1*state[7] + 3*b**2*k_off_1*state[9] + 2*Ca*k_on_2*state[2] + 2*b*k_off_2*state[14])

        def dMature_31(state, Ca):

            return(-1*(2*Ca*k_on_1 + 3*b**2*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**3*s + k_unmat)*state[9] + 3*Ca*k_on_1*state[8] + 4*b**3*k_off_1*state[10] + 2*Ca*k_on_2*state[3] + 2*b*k_off_2*state[15])

        def dMature_41(state, Ca):

            return(-1*(Ca*k_on_1 + 4*b**3*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**4*s + k_unmat)*state[10] + 2*Ca*k_on_1*state[9] + 5*b**3*k_off_1*state[11] + 2*Ca*k_on_2*state[4] + 2*b*k_off_2*state[16])

        def dMature_51(state, Ca):

            return(-1*(5*b**4*k_off_1 + Ca*k_on_2 + k_off_2 + k_fuse_basal*f**5*s + k_unmat)*state[11] + Ca*k_on_1*state[10] + 2*Ca*k_on_2*state[5] + 2*b*k_off_2*state[17])

        ###

        def dMature_02(state, Ca):

            return(-1*(5*Ca*k_on_1 + 2*b*k_off_2 + k_fuse_basal*s**2 + k_unmat)*state[12] + k_off_1*state[13] + Ca*k_on_2*state[6])

        def dMature_12(state, Ca):

            return(-1*(4*Ca*k_on_1 + k_off_1 + 2*b*k_off_2 + k_fuse_basal*f*s**2 + k_unmat)*state[13] + 5*Ca*k_on_1*state[12] + 2*b*k_off_1*state[14] + Ca*k_on_2*state[7])

        def dMature_22(state, Ca):

            return(-1*(3*Ca*k_on_1 + 2*b*k_off_1 + 2*b*k_off_2 + k_fuse_basal*f**2*s**2 + k_unmat)*state[14] + 4*Ca*k_on_1*state[13] + 3*b**2*k_off_1*state[15] + Ca*k_on_2*state[8])

        def dMature_32(state, Ca):

            return(-1*(2*Ca*k_on_1 + 3*b**2*k_off_1 + 2*b*k_off_2 + Ca*k_on_2 + k_fuse_basal*f**3*s**2 + k_unmat)*state[15] + 3*Ca*k_on_1*state[14] + 4*b**3*k_off_1*state[16] + Ca*k_on_2*state[9])

        def dMature_42(state, Ca):

            return(-1*(Ca*k_on_1 + 4*b**3*k_off_1 + 2*b*k_off_2 + Ca*k_on_2 + k_fuse_basal*f**4*s**2 + k_unmat)*state[16] + 2*Ca*k_on_1*state[15] + 5*b**3*k_off_1*state[17] + Ca*k_on_2*state[10])

        def dMature_52(state, Ca):

            return(-1*(5*b**4*k_off_1 + 2*b*k_off_2 + k_fuse_basal*f**5*s**2 + k_unmat)*state[17] + Ca*k_on_1*state[16] + Ca*k_on_2*state[11])

        def dFused(state):

            return(k_fuse_basal*(state[0] + f*state[1] + f**2*state[2] + f**3*state[3] + f**4*state[4] + f**5*state[5] + s*state[6] + f*s*state[7] + f**2*s*state[8] + f**3*s*state[9] + f**4*s*state[10] + f**5*s*state[11] + s**2*state[12] + f*s**2*state[13] + f**2*s**2*state[14] + f**3*s**2*state[15] + f**4*s**2*state[16] + f**5*s**2*state[17]))

        ###

        def partials(t, y):
            state = y[0:-3]
            immature = y[-3]
            fused = y[-2]
            Ca = y[-1]

            return(dMature_00(state, Ca), dMature_10(state, Ca), dMature_20(state, Ca), dMature_30(state, Ca), dMature_40(state, Ca), dMature_50(state, Ca), dMature_01(state, Ca), dMature_11(state, Ca), dMature_21(state, Ca), dMature_31(state, Ca), dMature_41(state, Ca), dMature_51(state, Ca), dMature_02(state, Ca), dMature_12(state, Ca), dMature_22(state, Ca), dMature_32(state, Ca), dMature_42(state, Ca), dMature_52(state, Ca), dImmature(state), dFused(state), dCa(t, Ca))


        max_time = stimulus_times[-1] + stimulus_times[0]*2
        ts = np.linspace(0, max_time, int(max_time+1))

        self.sol = solve_ivp(partials, [0, max_time], state_ini, t_eval = ts)

        self.state_ini = state_ini
        self.ts = ts
        self.max_time = max_time
        self.n_sites = n_sites
        self.b = b
