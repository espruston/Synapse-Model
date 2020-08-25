import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from math import e

if __name__ == "__main__":


    K_D_3 = 3e5 #M^-1ms^-1 Hui
    K_D_7 = 7.333e3 #M^-1ms^-1 Brandt/Knight

    k_docking =
    k_undocking =
    k_maturation =
    k_dematuration =

    CDR =
    Facil =

    C_3 =
    C_7 =

    p_immature =
    p_mature =

    Ca_rest = 5e-8 #M
    Ca_spike = 2e-5 #M
    Ca_residual = 250e-9 #M
    T_Ca_decay = 40 #ms

    delta_t = 1e-2 #ms

    t_SS = 10000 #ms
    ts_SS = np.linspace(0, t_SS, t_SS*delta_t + 1)

    SS = [0 1 0] #start all vesicles in immature docked state
    for t in ts_SS:

    stimulus_times = [0 10]

    ts = np.linspace()
