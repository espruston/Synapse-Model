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

        syt17 = np.asarray(membrane_bound_syt1[1:])*np.asarray(membrane_bound_syt7[1:])
        #syt17 /= max(syt17[t0:t0+int(20/delta_t)]) #normalize to first peak
        Fused = [0]
        E_A = -40 #kT
        E_syt1 = 20
        E_syt7 = 2

        m1 = max(membrane_bound_syt1)
        m7 = max(membrane_bound_syt7)

        for i in range(len(ts)):

            k_fuse = e**(E_A + E_syt1*membrane_bound_syt1[i]/m1 + E_syt7*membrane_bound_syt7[i]/m7)
            Fused.append(Fused[-1] + k_fuse*delta_t)

        self.ts = ts
        self.Ca = Ca_sim
        self.syt1 = np.asarray(membrane_bound_syt1[1:])
        self.syt7 = np.asarray(membrane_bound_syt7[1:])
        self.syt1Ca = Ca_bound_syt1[1:]
        self.syt7Ca = Ca_bound_syt7[1:]
        self.Fused = Fused[1:]
        self.dFused = np.diff(Fused)


class three_sensor(object):

    def __init__(self, K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times):

        sigma = FWHM/2.35 #variance
        mu = 2*FWHM #time at which Ca_spike is maximal (ms

        ts = np.arange(-1000,max_time+delta_t, delta_t) #simulate from t = -1000 to t = max_time (inclusive) with a resolution of delta_t ms
        #simulate calcium influx as a gaussian peak with exponential decay of residual Ca
        Ca_sim = np.zeros(len(ts))

        Ca_sim[:] += Ca_rest
        t0 = int(-1*ts[0]/delta_t)
        for t in stimulus_times:

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
        self.sty17 = syt17
        self.syt1Ca = Ca_bound_syt1[1:]
        self.syt3Ca = Ca_bound_syt3[1:]
        self.syt7Ca = Ca_bound_syt7[1:]

class two_pool_three_sensor_1(object):
    #two pool model where syt3 acts like syt1 for the "distant" pool

    def __init__(self, size_1, size_2, k_1_basal, k_2_basal, K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times, **kwargs):

        E_A = kwargs.get('E_A', -40)
        E_syt1 = kwargs.get('E_syt1', 20)
        E_syt3 = kwargs.get('E_syt3', 20)
        E_syt7 = kwargs.get('E_syt7', 2)
        k_fuse_basal = kwargs.get('k_fuse_basal', 3.5e-7)
        f = kwargs.get('f', 28**5)
        s = kwargs.get('s', 510**2)

        sigma = FWHM/2.35 #variance
        mu = 2*FWHM #time at which Ca_spike is maximal (ms

        ts = np.arange(-1000,max_time+delta_t, delta_t) #simulate from t = -1000 to t = max_time (inclusive) with a resolution of delta_t ms
        #simulate calcium influx as a gaussian peak with exponential decay of residual Ca
        Ca_local = np.zeros(len(ts))
        Ca_res = np.zeros(len(ts))

        Ca_res[:] += Ca_rest

        for t in stimulus_times:
            t0 = int(-1*ts[0]/delta_t)
            spike_start_index = t0 + int(t/delta_t) #index for the Ca_sim vector corresponding to t
            spike_peak_index = t0 + int((t+mu)/delta_t) #index for the Ca_sim vector corresponding to the peak of the spike initiated at t

            Ca_local[spike_start_index:] += Ca_spike*np.exp(-1*((ts[t0:t0+len(ts)-spike_start_index] - mu)/sigma)**2/2)
            Ca_res[spike_start_index:] += Ca_residual*np.exp(-1*ts[t0:t0+len(ts)-spike_start_index]/T_Ca_decay) #when the pulse hits its peak, add in a residual decay component from that point on

        Ca_bound_syt1 = [0]
        Ca_bound_syt3 = [0]
        Ca_bound_syt7_1 = [0] #Ca bound syt7 in pool 1
        Ca_bound_syt7_2 = [0] #Ca bound syt7 in pool 2

        membrane_bound_syt1 = [0]
        membrane_bound_syt3 = [0]
        membrane_bound_syt7_1 = [0] #membrane bound syt7 in pool 1
        membrane_bound_syt7_2 = [0] #membrane bound syt7 in pool 2

        for i in range(len(ts)):

            Ca_1 = Ca_local[i]+Ca_res[i]
            Ca_2 = Ca_res[i]

            #for pool 1 (close pool)
            membrane_bound_syt1.append(membrane_bound_syt1[-1] + ((1 - membrane_bound_syt1[-1])*k_on_1*Ca_bound_syt1[-1] - membrane_bound_syt1[-1]*k_off_1)*delta_t)
            membrane_bound_syt7_1.append(membrane_bound_syt7_1[-1] + ((1 - membrane_bound_syt7_1[-1])*k_on_7*Ca_bound_syt7_1[-1] - membrane_bound_syt7_1[-1]*k_off_7)*delta_t)

            #for pool 2 (distant pool)
            membrane_bound_syt3.append(membrane_bound_syt3[-1] + ((1 - membrane_bound_syt3[-1])*k_on_3*Ca_bound_syt3[-1] - membrane_bound_syt3[-1]*k_off_3)*delta_t)
            membrane_bound_syt7_2.append(membrane_bound_syt7_2[-1] + ((1 - membrane_bound_syt7_2[-1])*k_on_7*Ca_bound_syt7_2[-1] - membrane_bound_syt7_2[-1]*k_off_7)*delta_t)


            Ca_bound_syt1.append(Ca_1**2/(K_D_1**2 + Ca_1**2))
            Ca_bound_syt3.append(Ca_2**2/(K_D_3**2 + Ca_2**2))
            Ca_bound_syt7_1.append(Ca_1**2/(K_D_7**2 + Ca_1**2))
            Ca_bound_syt7_2.append(Ca_2**2/(K_D_7**2 + Ca_2**2))

        membrane_bound_syt1_norm = np.asarray(membrane_bound_syt1)/max(membrane_bound_syt1)
        membrane_bound_syt3_norm = np.asarray(membrane_bound_syt3)/max(membrane_bound_syt3)
        membrane_bound_syt7_1_norm = np.asarray(membrane_bound_syt7_1)/max(max(membrane_bound_syt7_1),max(membrane_bound_syt7_2))
        membrane_bound_syt7_2_norm = np.asarray(membrane_bound_syt7_2)/max(max(membrane_bound_syt7_1),max(membrane_bound_syt7_2))

        #max_time ms simulation, time resolution of delta_t (ms)
        pool_1 = [size_1]
        pool_2 = [size_2]

        dFused_1 = []
        dFused_2 = []

        for i in range(len(ts)):

            # k_fuse_1 =
            # k_fuse_2 = e**(E_A + E_syt3*membrane_bound_syt3_norm[i] + E_syt7*membrane_bound_syt7_2_norm[i])

            k_17_1 =  membrane_bound_syt1[i]*membrane_bound_syt7_1[i]*e**(E_A + E_syt1 + E_syt7)
            k_1_1 = membrane_bound_syt1[i]*(1-membrane_bound_syt7_1[i])*e**(E_A + E_syt1)
            k_7_1 = membrane_bound_syt7_1[i]*(1-membrane_bound_syt1[i])*e**(E_A + E_syt7)

            dFused_1.append(pool_1[-1]*2e12*(k_17_1 + k_1_1 + k_7_1)*delta_t)
            pool_1.append(pool_1[-1] - dFused_1[-1] + (size_1 - pool_1[-1])*k_1_basal*delta_t)

            k_37_2 =  membrane_bound_syt3[i]*membrane_bound_syt7_2[i]*e**(E_A + E_syt3 + E_syt7)
            k_3_2 = membrane_bound_syt3[i]*(1-membrane_bound_syt7_2[i])*e**(E_A + E_syt3)
            k_7_2 = membrane_bound_syt7_2[i]*(1-membrane_bound_syt3[i])*e**(E_A + E_syt7)

            dFused_2.append(pool_2[-1]*2e12*(k_37_2 + k_3_2 + k_7_2)*delta_t)
            pool_2.append(pool_2[-1] - dFused_2[-1] + (size_2 - pool_2[-1])*k_2_basal*delta_t)


        self.ts = ts
        self.Ca_local = np.asarray(Ca_local)
        self.Ca_res = np.asarray(Ca_res)
        self.pool_1 = pool_1[1:]
        self.pool_2 = pool_2[1:]
        self.dFused_1 = np.asarray(dFused_1)
        self.dFused_2 = np.asarray(dFused_2)
        self.syt1 = np.asarray(membrane_bound_syt1[1:])
        self.syt3 = np.asarray(membrane_bound_syt3[1:])
        self.syt7_1 = np.asarray(membrane_bound_syt7_1[1:])
        self.syt7_2 = np.asarray(membrane_bound_syt7_2[1:])
        self.syt1Ca = Ca_bound_syt1[1:]
        self.syt3Ca = Ca_bound_syt3[1:]
        self.syt7_1Ca = Ca_bound_syt7_1[1:]
        self.syt7_2Ca = Ca_bound_syt7_2[1:]


class two_pool_three_sensor_2(object):
    #two pool model where syt3 acts as a refill enhancer

    def __init__(self, size_1, size_2, k_1_basal, k_2_basal, K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times, **kwargs):

        E_A = kwargs.get('E_A', -40)
        E_syt1 = kwargs.get('E_syt1', 20)
        E_syt3 = kwargs.get('E_syt3', 20)
        E_syt7 = kwargs.get('E_syt7', 2)
        k_fuse_basal = kwargs.get('k_fuse_basal', 3.5e-7)
        f = kwargs.get('f', 28**5)
        s = kwargs.get('s', 510**2)

        sigma = FWHM/2.35 #variance
        mu = 2*FWHM #time at which Ca_spike is maximal (ms

        ts = np.arange(-1000,max_time+delta_t, delta_t) #simulate from t = -1000 to t = max_time (inclusive) with a resolution of delta_t ms
        #simulate calcium influx as a gaussian peak with exponential decay of residual Ca
        Ca_local = np.zeros(len(ts))
        Ca_res = np.zeros(len(ts))

        Ca_res[:] += Ca_rest
        t0 = int(-1*ts[0]/delta_t)
        for t in stimulus_times:

            spike_start_index = t0 + int(t/delta_t) #index for the Ca_sim vector corresponding to t
            spike_peak_index = t0 + int((t+mu)/delta_t) #index for the Ca_sim vector corresponding to the peak of the spike initiated at t

            Ca_sim[spike_start_index:] += Ca_spike*np.exp(-1*((ts[t0:t0+len(ts)-spike_start_index] - mu)/sigma)**2/2)
            Ca_sim[spike_peak_index:] += Ca_residual*np.exp(-1*ts[t0:t0+len(ts)-spike_peak_index]/T_Ca_decay) #when the pulse hits its peak, add in a residual decay component from that point on

        Ca_bound_syt1 = [0]
        Ca_bound_syt3 = [0]
        Ca_bound_syt7 = [0]

        membrane_bound_syt1 = [0]
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


class two_pool_three_sensor_3(object):
    #two pool model where syt3 acts as a fusion enhancer for the distant pool

    def __init__(self, size_1, size_2, k_1_basal, k_2_basal, K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times, **kwargs):

        E_A = kwargs.get('E_A', -40)
        E_syt1 = kwargs.get('E_syt1', 20)
        E_syt3 = kwargs.get('E_syt3', 20)
        E_syt7 = kwargs.get('E_syt7', 2)
        k_fuse_basal = kwargs.get('k_fuse_basal', 3.5e-7)
        f = kwargs.get('f', 28**5)
        s = kwargs.get('s', 510**2)

        sigma = FWHM/2.35 #variance
        mu = 2*FWHM #time at which Ca_spike is maximal (ms

        ts = np.arange(-1000,max_time+delta_t, delta_t) #simulate from t = -1000 to t = max_time (inclusive) with a resolution of delta_t ms
        #simulate calcium influx as a gaussian peak with exponential decay of residual Ca
        Ca_local = np.zeros(len(ts))
        Ca_res = np.zeros(len(ts))

        Ca_res[:] += Ca_rest
        t0 = int(-1*ts[0]/delta_t)
        for t in stimulus_times:

            spike_start_index = t0 + int(t/delta_t) #index for the Ca_sim vector corresponding to t
            spike_peak_index = t0 + int((t+mu)/delta_t) #index for the Ca_sim vector corresponding to the peak of the spike initiated at t

            Ca_sim[spike_start_index:] += Ca_spike*np.exp(-1*((ts[t0:t0+len(ts)-spike_start_index] - mu)/sigma)**2/2)
            Ca_sim[spike_peak_index:] += Ca_residual*np.exp(-1*ts[t0:t0+len(ts)-spike_peak_index]/T_Ca_decay) #when the pulse hits its peak, add in a residual decay component from that point on
