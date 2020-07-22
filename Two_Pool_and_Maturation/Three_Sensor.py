import numpy as np
from math import e

class Skyler_dual_sensor(object):

    def __init__(self, K_D_1, K_D_7, k_on_1, k_on_7, k_off_1, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times):

        sigma = FWHM/2.35 #variance
        mu = 2*FWHM #time at which Ca_spike is maximal (ms)

        ts = np.arange(-1000,max_time+delta_t, delta_t) #simulate from t = -1000 to t = max_time (inclusive) with a resolution of delta_t ms
        #simulate calcium influx as a gaussian peak with exponential decay of residual Ca
        Ca_sim = np.zeros(len(ts))

        Ca_sim[:] += Ca_rest
        for t in stimulus_times:
            t0 = int(-1*ts[0]/delta_t)
            spike_start_index = t0 + int(t/delta_t) #index for the Ca_sim vector corresponding to t
            spike_peak_index = t0 + int((t+mu)/delta_t) #index for the Ca_sim vector corresponding to the peak of the spike initiated at t

            Ca_sim[spike_start_index:] += Ca_spike*np.exp(-1*((ts[t0:t0+len(ts)-spike_start_index] - mu)/sigma)**2/2)
            Ca_sim[spike_peak_index:] += Ca_residual*np.exp(-1*ts[t0:t0+len(ts)-spike_peak_index]/T_Ca_decay) #when the pulse hits its peak, add in a residual decay component from that point on

        Ca_bound_syt1 = [0]
        Ca_bound_syt7 = [0]

        membrane_bound_syt1 = [0]
        membrane_bound_syt7 = [0]

        #max_time+1000 ms simulation, time resolution of delta_t (ms) with steady state determination for Ca_rest
        for i in range(len(ts)):

            Ca = Ca_sim[i]

            membrane_bound_syt1.append(membrane_bound_syt1[-1] + ((1 - membrane_bound_syt1[-1])*k_on_1*Ca_bound_syt1[-1] - membrane_bound_syt1[-1]*k_off_1)*delta_t)#increases according to free_syt1*k_on*Ca_bound, decreases according to bound*k_off over delta_t
            membrane_bound_syt7.append(membrane_bound_syt7[-1] + ((1 - membrane_bound_syt7[-1])*k_on_7*Ca_bound_syt7[-1]- membrane_bound_syt7[-1]*k_off_7)*delta_t)

            Ca_bound_syt1.append(Ca**2/(K_D_1**2 + Ca**2))
            Ca_bound_syt7.append(Ca**2/(K_D_7**2 + Ca**2))




        self.ts = ts
        self.Ca = Ca_sim
        self.syt1 = membrane_bound_syt1[1:]
        self.syt7 = membrane_bound_syt7[1:]
        self.syt1Ca = Ca_bound_syt1[1:]
        self.syt7Ca = Ca_bound_syt7[1:]


class three_sensor(object):

    def __init__(self, K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times):

        sigma = FWHM/2.35 #variance
        mu = 2*FWHM #time at which Ca_spike is maximal (ms

        ts = np.arange(-1000,max_time+delta_t, delta_t) #simulate from t = -1000 to t = max_time (inclusive) with a resolution of delta_t ms
        #simulate calcium influx as a gaussian peak with exponential decay of residual Ca
        Ca_sim = np.zeros(len(ts))

        Ca_sim[:] += Ca_rest
        for t in stimulus_times:
            t0 = int(-1*ts[0]/delta_t)
            spike_start_index = t0 + int(t/delta_t) #index for the Ca_sim vector corresponding to t
            spike_peak_index = t0 + int((t+mu)/delta_t) #index for the Ca_sim vector corresponding to the peak of the spike initiated at t

            Ca_sim[spike_start_index:] += Ca_spike*np.exp(-1*((ts[t0:t0+len(ts)-spike_start_index] - mu)/sigma)**2/2)
            Ca_sim[spike_peak_index:] += Ca_residual*np.exp(-1*ts[t0:t0+len(ts)-spike_peak_index]/T_Ca_decay) #when the pulse hits its peak, add in a residual decay component from that point on

        Ca_bound_syt1 = [0]
        Ca_bound_syt3 = [0]
        Ca_bound_syt7 = [0]

        membrane_bound_syt1 = [0]
        membrane_bound_syt3 = [0]
        membrane_bound_syt7 = [0]

        #max_time ms simulation, time resolution of delta_t (ms)
        for i in range(len(ts)):

            Ca = Ca_sim[i]

            membrane_bound_syt1.append(membrane_bound_syt1[-1] + ((1 - membrane_bound_syt1[-1])*k_on_1*Ca_bound_syt1[-1] - membrane_bound_syt1[-1]*k_off_1)*delta_t)
            membrane_bound_syt3.append(membrane_bound_syt3[-1] + ((1 - membrane_bound_syt3[-1])*k_on_3*Ca_bound_syt3[-1] - membrane_bound_syt3[-1]*k_off_3)*delta_t)
            membrane_bound_syt7.append(membrane_bound_syt7[-1] + ((1 - membrane_bound_syt7[-1])*k_on_7*Ca_bound_syt7[-1] - membrane_bound_syt7[-1]*k_off_7)*delta_t)

            Ca_bound_syt1.append(Ca**2/(K_D_1**2 + Ca**2))
            Ca_bound_syt3.append(Ca**2/(K_D_3**2 + Ca**2))
            Ca_bound_syt7.append(Ca**2/(K_D_7**2 + Ca**2))


        self.ts = ts
        self.Ca = Ca_sim
        self.syt1 = membrane_bound_syt1[1:]
        self.syt3 = membrane_bound_syt3[1:]
        self.syt7 = membrane_bound_syt7[1:]
        self.syt1Ca = Ca_bound_syt1[1:]
        self.syt3Ca = Ca_bound_syt3[1:]
        self.syt7Ca = Ca_bound_syt7[1:]


class three_sensor_simple(object):

    def __init__(self, stimulus_times, p_base, delta_k_3, delta_p_7, delta_3, delta_7, k_off_3, k_off_7,):

        syt3 = 0
        syt7 = 0

        n_vesicles = 1

        ts = np.arange(0,max_time+delta_t, delta_t) #simulate from t = 0 to t = max_time (inclusive) with a resolution of delta_t ms
        Ca_sim = np.zeros(len(ts))

        Ca_sim[:] += Ca_rest
        for t in stimulus_times:
            spike_start_index = int(t/delta_t) #index for the Ca_sim vector corresponding to t
            spike_peak_index = int((t+mu)/delta_t) #index for the Ca_sim vector corresponding to the peak of the spike initiated at t

            Ca_sim[spike_start_index:] += Ca_spike*np.exp(-1*((ts[0:len(ts)-spike_start_index] - mu)/sigma)**2/2)
            Ca_sim[spike_peak_index:] += Ca_residual*np.exp(-1*ts[0:-1*spike_peak_index]/T_Ca_decay) #when the pulse hits its peak, add in a residual decay component from that point on

class three_sensor_ultra_simple(object):

    def __init__(self, stimulus_times, k_refill_basal, p_base, delta_k_3, delta_p_7, k_on_3, k_on_7, k_off_3, k_off_7, delta_t):

        syt3vec = []
        syt7vec = []
        syt3 = 0
        syt7 = 0

        k_refillvec = []

        n_vesicles = 1
        n_vesiclesvec = []

        EPSC = []

        ISIs = np.diff(stimulus_times)
        ISIs = np.append(ISIs, ISIs[0])

        for i in range(len(stimulus_times)):

            EPSC.append(n_vesicles*(p_base + syt7*delta_p_7))
            n_vesicles -= EPSC[-1]

            syt3 += (1-syt3)*k_on_3 #syt steps up as a result of the pulse by a constant percent of Ca-unbound syt
            syt7 += (1-syt7)*k_on_7

            ts = np.arange(0, ISIs[i]/delta_t, delta_t)
            for t in ts:

                k_refill = k_refill_basal + syt3*delta_k_3 #k_refill is greater than k_refill_basal by delta_k_3 when all syt3 are Ca bound

                n_vesicles += (1-n_vesicles)*k_refill*delta_t

                syt3 *= 1-k_off_3*delta_t #Ca bound syt decays with rate constant over delta_t between stimuli
                syt7 *= 1-k_off_7*delta_t

                #if t%(delta_t*10) == 0:
                k_refillvec.append(k_refill)
                n_vesiclesvec.append(n_vesicles)
                syt3vec.append(syt3)
                syt7vec.append(syt7)

        self.syt3 = syt3vec
        self.syt7 = syt7vec
        self.k_refill = k_refillvec
        self.vesicles = n_vesiclesvec
        if sum(ISIs)/delta_t < 1e5:
            self.ts = np.arange(0, sum(ISIs)/delta_t, delta_t)*delta_t
        self.EPSC = np.asarray(EPSC)
        self.norm_val = EPSC[0]
