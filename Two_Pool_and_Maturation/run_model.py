import numpy as np
import pandas as pd
import os
from TwoPoolAndMat import *
from Three_Sensor import *
from matplotlib import pyplot as plt
from math import e

if __name__ == "__main__":

    data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', nrows = 20)
    data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 22, nrows = 20)
    data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 44, nrows = 20)
    data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', usecols = 'C,D', skiprows = 66, nrows = 20)

    EPSC_1hz = data_1hz.to_numpy()[:,0]
    stdev_1hz = data_1hz.to_numpy()[:,1]
    EPSC_10hz = data_10hz.to_numpy()[:,0]
    stdev_10hz = data_10hz.to_numpy()[:,1]
    EPSC_20hz = data_20hz.to_numpy()[:,0]
    stdev_20hz = data_20hz.to_numpy()[:,1]
    EPSC_50hz = data_50hz.to_numpy()[:,0]
    stdev_50hz = data_50hz.to_numpy()[:,1]

    #two_pool(r, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F)
    #maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)
    #maturation_2(stimulus_times, .002, 8e8, 120, 1.4e8, 4e3, 4e7, 2e3, 28, 510, 3.5e-4, 5e-8, 2.5e-5, 0.04)
    #Skyler_dual_sensor(K_D_1, K_D_7, k_on_1, k_on_7, k_off_1, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times)
    #Evan_dual_sensor(k_on_1, k_on_7, k_off_1, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times)
    #three_sensor(K_D_1, K_D_3, K_D_7, k_on_1, k_on_3, k_on_7, k_off_1, k_off_3, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times)
    #
    # r = 10
    # n_pulses = 20
    #
    '''two pool parameters'''
    # size_fast = .68
    # size_slow = 1-size_fast
    # p_fast = .02
    # p_slow = .26
    # T_fast = .11
    # T_slow = 4.9
    # delta_F = 0
    # T_F = 0.00001
    #
    '''maturation parameters'''
    # p_immature = .2
    # p_mature = 1-.198
    # p_facilitated = 1
    # T_refill = 12
    # T_maturation = 250
    # T_facilitation = 1e30
    #
    '''dual&three sensor parameters'''
    #K_Ds for Ca binding of syts
    # K_D_1 =
    K_D_1 = 41 #syt1 K_D for membrane binding, uM, Brandt/Knight
    K_D_3 = 2 #syt3, uM, Sugita
    K_D_7 = 1.5 #syt 7, uM, Knight

    K_A_1 = 43e-6 #half max Ca concentration for membrane binding of syt1, M, Knight
    K_A_3 = 0
    K_A_7 = 5.4e-6 #syt7, M, Knight
    #
    #membrane association constants
    #Ca dependent k_on
    #k_on_1 = 1.63e5 #syt1 M-1ms-1, Knight
    #k_on_1 = 1.2e7 #Hui
    #k_on_3 = 3e5 #syt3 M-1ms-1, Hui
    #k_on_7 = 7.333e3 #syt7 M-1ms-1, Knight
    #k_on_7 = 2.333e4 #fits skylers model well with syt1 from Knight
    #k_on_7 = 3e5 #estimate from likeness to syt3 Hui
    #
    # #Ca indpendent k_on
    k_on_1 = .340 #ms-1, Jackman
    # k_on_3 = .190 #ms-1, est. from syt 7 vals, Jackman
    k_on_7 = .157 #ms-1, Jackman
    #
    k_off_1 = .670 #syt1 ms-1, 90% amplitude, major rate, double exponential, Knight
    # #k_off_1 = .378 #Hui
    # k_off_3 = .0128 #syt3 ms-1, Hui
    k_off_7 = .011 #syt7 ms-1, Knight
    # #k_off_7 = .019*.25 + .008*.75 #Hui, crude approx of the double exponential
    #
    delta_t = 1e-1 #time resolution, ms
    max_time = 500 #ms
    #stimulus_times = np.arange(0,200,20)
    stimulus_times = [0,20]
    Ca_rest = 5e-2 #Resting calcium uM, Jackman
    Ca_residual = 25e-2 #Residual calcium uM, Jackman
    T_Ca_decay = 40 #Residual calcium decay constant ms, Jackman
    Ca_spike = 25 #Local calcium after pulse uM, Jackman
    FWHM = .34 #Local calcium full width half maximum ms

    output = Skyler_dual_sensor(K_D_1, K_D_7, k_on_1, k_on_7, k_off_1, k_off_7, Ca_rest, Ca_residual, T_Ca_decay, Ca_spike, FWHM, delta_t, max_time, stimulus_times)

    #stimulus_times = np.arange(0,100,10)

    #output = three_sensor_ultra_simple(stimulus_times, k_refill_basal, p_base, delta_k_3, delta_p_7, k_on_3, k_on_7, k_off_3, k_off_7, delta_t)
    #output = three_sensor_ultra_simple(stimulus_times, 3.5e-6, .8, 5e-7, 0, 0, 0, .00001, .0019, 1e-2)
    #print(output.EPSC[-1])
    #plt.plot(stimulus_times, output.EPSC/output.norm_val, 'bo')
    #plt.plot(output.ts, output.syt3, label = "syt3")
    #plt.plot(output.ts, output.syt7, label = "syt7")
    #plt.plot(output.ts, output.vesicles)
    # plt.ylabel("ESPSC norm")
    # plt.xlabel("Time (ms)")
    # plt.xlim(0,stimulus_times[-1]+100)
    # plt.ylim(0,1.2)
    #plt.legend()

    plt.plot(output.ts, output.Ca)
    plt.ylabel("Ca concentration (uM)")
    plt.xlabel("time (ms)")
    plt.title("Simulated presynaptic Ca for 50hz PPR")

    plt.xlim(-10,max_time)
    plt.yscale('log')

    #plt.show()
    #
    # plt.plot(output.ts, output.syt1/output.norm_val, label = "Membrane bound Syt 1")
    # plt.plot(output.ts, output.syt7/output.norm_val, label = "Membrane bound Syt 7")
    # plt.ylabel("Bound isoform (norm.)")
    # plt.xlabel("time (ms)")
    # plt.title("Simulated SYT membrane binding (single pulse)")
    # plt.xlim(0,max_time)
    # plt.ylim(0,1)

    #For Skyler_dual_sensor
    #PPR plot
    # t0 = int(1000/delta_t)
    # syt17 = np.asarray(output.syt1[1:])*np.asarray(output.syt7[1:])
    # syt17 /= max(syt17[t0:t0+int(20/delta_t)])
    # plt.plot(output.ts, syt17)
    # plt.ylabel("Syt1*Syt7 (norm)")
    # plt.xlabel("time (ms)")
    # plt.title("Simulated multiplicative membrane binding PPR")
    # plt.xlim(-10,max_time)

    #plt.plot(output.ts, output.syt1[1:]/max(output.syt1), label = "Membrane bound Syt 1")
    #plt.plot(output.ts, output.syt7[1:]/max(output.syt7), label = "Membrane bound Syt 7")
    #plt.plot(output.ts, output.syt1[1:], label = "Membrane bound Syt 1") #only normalize syt1 ???
    #plt.plot(output.ts, output.syt7[1:], label = "Membrane bound Syt 7")
    #plt.plot(output.ts, output.syt1Ca[1:]/output.norm_val_Ca, label = "Ca bound Syt 1")
    #plt.plot(output.ts, output.syt7Ca[1:]/output.norm_val_Ca, label = "Ca bound Syt 7")
    # plt.ylabel("Bound isoform (norm.)")
    # plt.xlabel("time (ms)")
    # plt.title("Simulated SYT membrane binding (single pulse)")
    # plt.xlim(-10,max_time)

    #plt.legend()
    #
    plt.show()

    '''usage'''

    #stimulus_times = np.linspace(20,200,10)
    #output = maturation_2(stimulus_times, .002, 8e8, 120, 1.4e8, 4e3, 4e7, 2e3, 28, 510, 3.5e-4, 5e-8, 2.5e-5, 0.04)

    #plt.plot(output.ts, output.sol.y)
    #plt.show()

    # # output1 = two_pool(r, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F)
    # #
    # output2 = maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)
    #
    # # fig = plt.plot(range(n_pulses), output1.EPSC, label = "two pool")
    # # plt.plot(range(n_pulses), output2.EPSC, label = "maturation")
    # # plt.scatter(range(n_pulses), EPSC_10hz, label = "data")
    # # plt.ylim(0, max([max(output1.EPSC), max(output2.EPSC)]))
    # # plt.xlim(0,n_pulses)
    #
    # fig, ax = plt.subplots()
    # plt.subplots_adjust(left=0.25, bottom=0.25)
    #
    # l, = plt.plot(range(n_pulses), output2.EPSC, label = "EPSC")
    # m, = plt.plot(range(n_pulses), output2.EPSC_immature, label = "Immature")
    # n, = plt.plot(range(n_pulses), output2.EPSC_mature, label = "Mature")
    # o, = plt.plot(range(n_pulses), output2.EPSC_facilitated, label = "Facilitated")
    # plt.xlim(0,n_pulses)
    # plt.ylim(0,1)
    # plt.scatter(range(n_pulses), EPSC_20hz, label = "data")
    # fig.tight_layout()
    # plt.legend()
    # plt.show()
