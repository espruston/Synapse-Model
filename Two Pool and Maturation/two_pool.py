import numpy as np
from matplotlib import pyplot as plt
from math import e

if __name__ == "__main__":

    r = 10

    m = .68
    fastpool_initial = m
    slowpool_initial = (1-m)

    fastpool = fastpool_initial
    slowpool = slowpool_initial

    p_fast = 0.02
    p_slow = .26

    T_fast = .11
    T_slow = 4.9

    F = 1
    K_F = 1.5
    delta_F = .5
    T_F = 0.07
    BG_F = 1e-5
    sat_F = K_F
    facil = BG_F

    n_pulses = 20

    EPSC = []
    fastEPSC = []
    slowEPSC = []

    # for i in np.arange(n_pulses):
    #
    #     facil = facil*e**(-1/(T_F*r))
    #     F = 1 + sat_F*facil/(K_F + facil)
    #
    #     slowrecover = (slowpool_initial - slowpool)*(1-e**(-1/(r*T_slow)))
    #     slowpool += slowrecover
    #     slowEPSC.append(slowpool*p_slow*F)
    #     slowpool -= slowpool*p_slow*F
    #
    #     fastrecover = (fastpool_initial - fastpool)*(1-e**(-1/(r*T_fast)))
    #     fastpool += fastrecover
    #     fastEPSC.append(fastpool*p_fast*F)
    #     fastpool -= fastpool*p_fast*F
    #
    #     EPSC.append(slowEPSC[-1] + fastEPSC[-1])
    #
    #     facil += delta_F

    for i in range(n_pulses):

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
        facil = facil*e**(-1/(T_F*r))
        F = 1 + sat_F*facil/(K_F + facil)

    fig = plt.plot(np.arange(n_pulses), np.asarray(EPSC)/EPSC[0], label = "EPSC")
    plt.plot(np.arange(n_pulses), np.asarray(fastEPSC)/EPSC[0], label = "Fast pool")
    plt.plot(np.arange(n_pulses), np.asarray(slowEPSC)/EPSC[0], label = "Slow pool")
    plt.xlim(0,n_pulses)
    plt.legend()
    plt.show()
