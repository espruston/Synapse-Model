import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
if __name__ == "__main__":

    K_D_1 = 43e-6 #syt1 Ca K_D, M, Brandt/Knight
    K_D_3 = 2e-6 #syt3 Ca K_D, M, Sugita
    K_D_7 = 1.5e-6 #syt7 Ca K_D, M, Knight

    #membrane association constants
    #Ca dependent k_on
    #k_on_1 = 1.63e5 #syt1 M-1ms-1, Knight
    #k_on_1 = 1.2e7 #Hui
    #k_on_3 = 3e5 #syt3 M-1ms-1, Hui
    #k_on_7 = 7.333e3 #syt7 M-1ms-1, Knight
    #k_on_7 = 2.333e4 #fits skylers model well with syt1 from Knight
    #k_on_7 = 3e5 #estimate from likeness to syt3 Hui

    #Ca indpendent k_on
    k_on_1 = .415 #ms-1, Jackman
    k_on_3 = .190 #ms-1, est. from syt 7 vals, Jackman
    k_on_7 = .190 #ms-1, Jackman

    k_off_1 = .670 #syt1 ms-1, 90% amplitude, major rate, double exponential, Knight
    #k_off_1 = .378 #Hui
    k_off_3 = .0128 #syt3 ms-1, Hui
    k_off_7 = .011 #syt7 ms-1, Knight
    #k_off_7 = .019*.25 + .008*.75 #Hui, crude approx of the double exponential

    cooperativity = 2

    Ca_rest = 50e-9 #Resting calcium M, Jackman
    Ca_residual = 250e-9 #Residual calcium M, Jackman
    T_Ca_decay = 40 #Residual calcium decay constant s, Jackman
    Ca_spike = 25e-6 #Local calcium after pulse M, Jackman
    FWHM = .34 #Local calcium full width half maximum ms
    sigma = FWHM/2.35 #variance
    mu = 2*FWHM #time at which Ca_spike is maximal (ms)

    delta_t = 1e-3 #time resolution, ms

    def SinglePulse(sigma, mu, T_Ca_decay, Ca_rest):

        t_spike = np.arange(0, mu+3*sigma+delta_t, delta_t)
        spike_sim = Ca_spike*np.exp(-1*((t_spike - mu)/sigma)**2/2) + Ca_rest

        t_residual = np.arange(delta_t, 500-t_spike[-1], delta_t)
        residual_sim = (spike_sim[-1] - Ca_rest)*np.exp(-1*t_residual/T_Ca_decay) + Ca_rest

        Ca_sim = np.concatenate((spike_sim, residual_sim))
        ts = np.concatenate((t_spike, t_residual+t_spike[-1]))

        return Ca_sim, ts

    def PPR(sigma, mu, T_Ca_decay, Ca_rest):

        t_spike = np.arange(0, mu+3*sigma+delta_t, delta_t) #simulate the spike up until after 3 sigma
        spike_sim = Ca_spike*np.exp(-1*((t_spike - mu)/sigma)**2/2) + Ca_rest

        t_residual1 = np.arange(delta_t, 20, delta_t)
        residual_sim1 = (spike_sim[-1] - Ca_rest)*np.exp(-1*t_residual1/T_Ca_decay) + Ca_rest
        residual_sim1[-1*len(spike_sim):] += spike_sim

        t_residual2 = np.arange(delta_t, 500-t_residual1[-1]-t_spike[-1], delta_t)
        residual_sim2 = (residual_sim1[-1] - Ca_rest)*np.exp(-1*t_residual2/T_Ca_decay) + Ca_rest

        Ca_sim = np.concatenate((spike_sim, residual_sim1, residual_sim2))
        ts = np.concatenate((t_spike, t_residual1+t_spike[-1], t_residual2+t_residual1[-1]+t_spike[-1]))

        return Ca_sim, ts

    Ca_sim, ts = SinglePulse(sigma, mu, T_Ca_decay, Ca_rest)

    # plt.plot(ts, Ca_sim)
    # plt.ylabel("Ca concentration (M)")
    # plt.xlabel("time (ms)")
    # plt.title("Simulated presynaptic Ca for 50hz PPR")
    # plt.xlim(-10,500)
    # plt.yscale('log')
    #
    # plt.show()
    #
    Ca_bound_syt1 = [0]
    Ca_bound_syt3 = [0]
    Ca_bound_syt7 = [0]

    bound_syt1 = [0]
    bound_syt3 = [0]
    bound_syt7 = [0]

    #500 ms simulation, time resolution of delta_t (ms)
    for i in range(len(ts)):

        Ca = Ca_sim[i]

        bound_syt1.append(bound_syt1[-1] + ((1-bound_syt1[-1])*k_on_1*Ca_bound_syt1[-1] - bound_syt1[-1]*k_off_1)*delta_t) #multiplied by delta_t to account for rates being in ms and timestep being in us
        bound_syt3.append(bound_syt3[-1] + ((1-bound_syt3[-1])*k_on_3*Ca_bound_syt3[-1] - bound_syt3[-1]*k_off_3)*delta_t)
        bound_syt7.append(bound_syt7[-1] + ((1-bound_syt7[-1])*k_on_7*Ca_bound_syt7[-1] - bound_syt7[-1]*k_off_7)*delta_t)

        Ca_bound_syt1.append(Ca**2/(K_D_1+Ca**2))
        Ca_bound_syt3.append(Ca**2/(K_D_3+Ca**2))
        Ca_bound_syt7.append(Ca**2/(K_D_7+Ca**2))

    norm_val = max(max(bound_syt1), max(bound_syt3), max(bound_syt7))
    #norm_val = max(max(bound_syt1), max(bound_syt7))
    norm_val_Ca = max(max(Ca_bound_syt1), max(Ca_bound_syt3), max(Ca_bound_syt7))

    # plt.plot(ts, bound_syt1[1:]/norm_val, label = "Syt 1")
    # plt.plot(ts, bound_syt3[1:]/norm_val, label = "Syt 3")
    # plt.plot(ts, bound_syt7[1:]/norm_val, label = "Syt 7")
    # plt.ylabel("Bound isoform (norm.)")
    # plt.xlabel("time (ms)")
    # plt.title("Simulated SYT membrane binding (single pulse)")
    # plt.xlim(0,500)
    # plt.ylim(0,1)
    # plt.legend()

    # plt.plot(ts, Ca_bound_syt1[1:]/norm_val_Ca, label = "Syt 1")
    # plt.plot(ts, Ca_bound_syt3[1:]/norm_val_Ca, label = "Syt 3")
    # plt.plot(ts, Ca_bound_syt7[1:]/norm_val_Ca, label = "Syt 7")
    # plt.ylabel("Bound isoform (norm.)")
    # plt.xlabel("time (ms)")
    # plt.title("Simulated SYT membrane binding (single pulse)")
    # plt.xlim(0,500)
    # plt.legend()

    plt.show()
