import pandas as pd
import numpy as np
from math import log, e
from matplotlib import pyplot as plt
from sklearn import metrics


if __name__ == "__main__":

    #import raw data to pandas dataframes
    CF_data_1hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 1, nrows = 20)
    CF_data_10hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 23, nrows = 20)
    CF_data_20hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 45, nrows = 20)
    CF_data_50hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 67, nrows = 20)

    PF_data_1hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 1, nrows = 20)
    PF_data_10hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 47, nrows = 20)
    PF_data_20hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 70, nrows = 20)
    PF_data_50hz = pd.read_excel(r'~\Dropbox\Work\Jackman Lab\Modeling\Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 93, nrows = 20)

    #define data series as numpy vectors
    CF_WT_1hz_EPSCs = CF_data_1hz.to_numpy()[:,0]
    CF_WT_1hz_SEMs = CF_data_1hz.to_numpy()[:,1]
    CF_KO_1hz_EPSCs = CF_data_1hz.to_numpy()[:,2]
    CF_KO_1hz_SEMs = CF_data_1hz.to_numpy()[:,3]
    CF_WTCaBlock_1hz_EPSCs = CF_data_1hz.to_numpy()[:,4]
    CF_WTCaBlock_1hz_SEMs = CF_data_1hz.to_numpy()[:,5]
    CF_KOCaBlock_1hz_EPSCs = CF_data_1hz.to_numpy()[:,6]
    CF_KOCaBlock_1hz_SEMs = CF_data_1hz.to_numpy()[:,7]

    CF_WT_10hz_EPSCs = CF_data_10hz.to_numpy()[:,0]
    CF_WT_10hz_SEMs = CF_data_10hz.to_numpy()[:,1]
    CF_KO_10hz_EPSCs = CF_data_10hz.to_numpy()[:,2]
    CF_KO_10hz_SEMs = CF_data_10hz.to_numpy()[:,3]
    CF_WTCaBlock_10hz_EPSCs = CF_data_10hz.to_numpy()[:,4]
    CF_WTCaBlock_10hz_SEMs = CF_data_10hz.to_numpy()[:,5]
    CF_KOCaBlock_10hz_EPSCs = CF_data_10hz.to_numpy()[:,6]
    CF_KOCaBlock_10hz_SEMs = CF_data_10hz.to_numpy()[:,7]

    CF_WT_20hz_EPSCs = CF_data_20hz.to_numpy()[:,0]
    CF_WT_20hz_SEMs = CF_data_20hz.to_numpy()[:,1]
    CF_KO_20hz_EPSCs = CF_data_20hz.to_numpy()[:,2]
    CF_KO_20hz_SEMs = CF_data_20hz.to_numpy()[:,3]
    CF_WTCaBlock_20hz_EPSCs = CF_data_20hz.to_numpy()[:,4]
    CF_WTCaBlock_20hz_SEMs = CF_data_20hz.to_numpy()[:,5]
    CF_KOCaBlock_20hz_EPSCs = CF_data_20hz.to_numpy()[:,6]
    CF_KOCaBlock_20hz_SEMs = CF_data_20hz.to_numpy()[:,7]

    CF_WT_50hz_EPSCs = CF_data_50hz.to_numpy()[:,0]
    CF_WT_50hz_SEMs = CF_data_50hz.to_numpy()[:,1]
    CF_KO_50hz_EPSCs = CF_data_50hz.to_numpy()[:,2]
    CF_KO_50hz_SEMs = CF_data_50hz.to_numpy()[:,3]
    CF_WTCaBlock_50hz_EPSCs = CF_data_50hz.to_numpy()[:,4]
    CF_WTCaBlock_50hz_SEMs = CF_data_50hz.to_numpy()[:,5]
    CF_KOCaBlock_50hz_EPSCs = CF_data_50hz.to_numpy()[:,6]
    CF_KOCaBlock_50hz_SEMs = CF_data_50hz.to_numpy()[:,7]

    PF_WT_1hz_EPSCs = PF_data_1hz.to_numpy()[:,0]
    PF_WT_1hz_SEMs = PF_data_1hz.to_numpy()[:,1]
    PF_KO_1hz_EPSCs = PF_data_1hz.to_numpy()[:,2]
    PF_KO_1hz_SEMs = PF_data_1hz.to_numpy()[:,3]
    PF_WTCaBlock_1hz_EPSCs = PF_data_1hz.to_numpy()[:,4]
    PF_WTCaBlock_1hz_SEMs = PF_data_1hz.to_numpy()[:,5]
    PF_KOCaBlock_1hz_EPSCs = PF_data_1hz.to_numpy()[:,6]
    PF_KOCaBlock_1hz_SEMs = PF_data_1hz.to_numpy()[:,7]

    PF_WT_10hz_EPSCs = PF_data_10hz.to_numpy()[:,0]
    PF_WT_10hz_SEMs = PF_data_10hz.to_numpy()[:,1]
    PF_KO_10hz_EPSCs = PF_data_10hz.to_numpy()[:,2]
    PF_KO_10hz_SEMs = PF_data_10hz.to_numpy()[:,3]
    PF_WTCaBlock_10hz_EPSCs = PF_data_10hz.to_numpy()[:,4]
    PF_WTCaBlock_10hz_SEMs = PF_data_10hz.to_numpy()[:,5]
    PF_KOCaBlock_10hz_EPSCs = PF_data_10hz.to_numpy()[:,6]
    PF_KOCaBlock_10hz_SEMs = PF_data_10hz.to_numpy()[:,7]

    PF_WT_20hz_EPSCs = PF_data_20hz.to_numpy()[:,0]
    PF_WT_20hz_SEMs = PF_data_20hz.to_numpy()[:,1]
    PF_KO_20hz_EPSCs = PF_data_20hz.to_numpy()[:,2]
    PF_KO_20hz_SEMs = PF_data_20hz.to_numpy()[:,3]
    PF_WTCaBlock_20hz_EPSCs = PF_data_20hz.to_numpy()[:,4]
    PF_WTCaBlock_20hz_SEMs = PF_data_20hz.to_numpy()[:,5]
    PF_KOCaBlock_20hz_EPSCs = PF_data_20hz.to_numpy()[:,6]
    PF_KOCaBlock_20hz_SEMs = PF_data_20hz.to_numpy()[:,7]

    PF_WT_50hz_EPSCs = PF_data_50hz.to_numpy()[:,0]
    PF_WT_50hz_SEMs = PF_data_50hz.to_numpy()[:,1]
    PF_KO_50hz_EPSCs = PF_data_50hz.to_numpy()[:,2]
    PF_KO_50hz_SEMs = PF_data_50hz.to_numpy()[:,3]
    PF_WTCaBlock_50hz_EPSCs = PF_data_50hz.to_numpy()[:,4]
    PF_WTCaBlock_50hz_SEMs = PF_data_50hz.to_numpy()[:,5]
    PF_KOCaBlock_50hz_EPSCs = PF_data_50hz.to_numpy()[:,6]
    PF_KOCaBlock_50hz_SEMs = PF_data_50hz.to_numpy()[:,7]

    ### CHANGE THESE TWO VALUES ###
    data = CF_WT_20hz_EPSCs
    hz = 20
    ### ----------------------- ###

    x = np.linspace(0,19*1000/hz,20)

    lny = []
    for i in range(len(data)):
        lny.append(log(data[i]))

    b = sum(x*lny)/sum(x**2)

    print('R squared:', metrics.r2_score(data, e**(b*x)))
    print('Equation: y = e^(',b,'* x)')

    fig, ax1 = plt.subplots()
    ax1.scatter((x+1000/hz)*hz/1000, data)
    ax1.set_xlabel('Pulse #')
    ax1.set_xlim(1,22)
    ax1.set_xticks(np.linspace(1,20,20))
    ax1.set_ylabel('$EPSC/EPSC_{1}$')
    ax1.set_ylim(0,max(data)+.5)

    xlong = np.linspace(0,19*1000/hz+2*1000/hz,22)
    ax2 = ax1.twiny()
    ax2.plot(xlong,e**(b*xlong))
    ax2.set_xlabel('t (ms)')
    ax2.set_xlim(0,max(xlong))
    ax2.set_xticks(xlong)

    plt.grid(axis = 'x')

    plt.show()
