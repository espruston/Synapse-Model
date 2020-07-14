import numpy as np
import pandas as pd
import os
from TwoPoolAndMat import two_pool, maturation, maturation_2
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

    stimulus_times = np.linspace(20,200,10)
    output = maturation_2(stimulus_times, .002, 8e8, 120, 1.4e8, 4e3, 4e7, 2e3, 28, 510, 3.5e-4, 5e-8, 2.5e-5, 0.04)

    plt.plot(output.ts, output.sol.y)
    plt.show()
    #two_pool(r, n_pulses, size_fast, size_slow, p_fast, p_slow, T_fast, T_slow, delta_F, T_F)
    #maturation(r, n_pulses, p_immature, p_mature, p_facilitated, T_refill, T_maturation, T_facilitation)
    #
    # r = 10
    # n_pulses = 20
    #
    # #two pool parameters
    # size_fast = .68
    # size_slow = 1-size_fast
    # p_fast = .02
    # p_slow = .26
    # T_fast = .11
    # T_slow = 4.9
    # delta_F = 0
    # T_F = 0.00001
    #
    # #maturation parameters
    # p_immature = .2
    # p_mature = 1-.198
    # p_facilitated = 1
    # T_refill = 12
    # T_maturation = 250
    # T_facilitation = 1e30
    #
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
