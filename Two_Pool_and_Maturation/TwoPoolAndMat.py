import numpy as np
from scipy.integrate import solve_ivp
from math import e, floor

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

# class sequentialRK1(object):
#
#     def __init__(self, stimulus_times, step_size, k_basal, k_prim, k_unprim, k_fill, k_unfill, k_on, k_off, b, gamma, T_Ca_decay):
#
#         ts = np.linspace(1,stimulus_times[-1],stimulus_times[-1])
#
#         Ca_1 = Ca_influx*(e**(-1*t/T_Ca_decay) #simnple exponential decay of calcium from max value of Ca_influx
#         Ca = np.zeros(stimulus_times[-1]) #functional alpha function
#
#         for i in range(len(stimulus_times)):
#             stimulus = stimulus_times[i]
#             Ca[stimulus-1:stimulus_times[-1]-1] += alpha_1[0:stimulus_times[-1]-stimulus]) #calculate effect of each stimulus on the Ca level
#
#         def vectorfield(states, t):
#
#             Unprime, Prime, Ca0_Releaseable, Ca1_Releaseable, Ca2_Releaseable, Ca3_Releaseable, Ca4_Releaseable, Ca5_Releaseable = states
#
#             if t in ts:
#                 Ca_val = Ca[t]
#             else:
#                 Ca_val = Ca[floor(t)]*e**(-1*(int(t) - floor(t))/T_Ca_decay)
#
#             #create f = [dUnprime, dPrime, dCa0_Releaseable, dCa1_Releaseable, dCa2_Releaseable, dCa3_Releaseable, dCa4_Releaseable, dCa5_Releaseable, V_fused]
#             f = [k_basal - k_prim*Ca_val*Unprime + k_unprim*Prime,
#             k_prim*Ca_val*Unprime - (k_unprim + k_fill)*Prime + k_unfill*Ca0_Releaseable,
#             k_fill*Prime - (5*k_on*Ca_val + k_unfill)*Ca0_Releaseable + k_off*b**0*Ca1_Releaseable, 5*k_on*Ca_val*Ca0_Releaseable - (1*k_off*b**0 + 4*k_on*Ca_val)*Ca1_Releaseable + 2*k_off*b**1*Ca2_Releaseable,
#             4*k_on*Ca_val*Ca1_Releaseable - (2*k_off*b**1 + 3*k_on*Ca_val)*Ca2_Releaseable + 3*k_off*b**2*Ca3_Releaseable,
#             3*k_on*Ca_val*Ca2_Releaseable - (3*k_off*b**2 + 2*k_on*Ca_val)*Ca3_Releaseable + 4*k_off*b**3*Ca4_Releaseable,
#             2*k_on*Ca_val*Ca3_Releaseable - (4*k_off*b**3 + 1*k_on*Ca_val)*Ca4_Releaseable + 5*k_off*b**4*Ca5_Releaseable,
#             k_on*Ca_val*Ca4_Releaseable - (gamma + 5*k_off*b**4)*Ca5_Releaseable,
#
#
#             return f
#
#         initial_conditions = [1, 1, .2, 0, 0, 0, 0, 0]
#         sols = solve_ivp(vectorfield, (1,t[-1]+100), initial_conditions, t_eval = ts)
