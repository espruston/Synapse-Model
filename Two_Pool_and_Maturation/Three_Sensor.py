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

class two_pool_three_sensor(object):

    def __init__(self, size_1, size_2, p_1, p_2, k_1_basal, k_2_basal, K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times):

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
