import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

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

    x = np.linspace(1,len(CF_WT_1hz_EPSCs),len(CF_WT_1hz_EPSCs))

    fig, axs = plt.subplots(4, 2, sharex = True, sharey = 'col')

    axs[0,0].errorbar(x, CF_WT_1hz_EPSCs, yerr = CF_WT_1hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[0,0].errorbar(x, CF_KO_1hz_EPSCs, yerr = CF_KO_1hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[0,0].errorbar(x, CF_WTCaBlock_1hz_EPSCs, yerr = CF_WTCaBlock_1hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[0,0].errorbar(x, CF_KOCaBlock_1hz_EPSCs, yerr = CF_KOCaBlock_1hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    axs[0,0].set_xlim(1, 21)
    plt.xticks(np.arange(min(x), max(x)+1, 1.0))
    axs[0,0].set_ylim(0,1.25)
    axs[0,0].set_title('Normalized CF EPSCs at 1 Hz')
    #axs[0].set_xlabel('Pulse #')
    axs[0,0].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    axs[0,0].legend(loc = "lower left", ncol = 2)

    axs[1,0].errorbar(x, CF_WT_10hz_EPSCs, yerr = CF_WT_10hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[1,0].errorbar(x, CF_KO_10hz_EPSCs, yerr = CF_KO_10hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[1,0].errorbar(x, CF_WTCaBlock_10hz_EPSCs, yerr = CF_WTCaBlock_10hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[1,0].errorbar(x, CF_KOCaBlock_10hz_EPSCs, yerr = CF_KOCaBlock_10hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    # axs[1].set_xlim(1, 21)
    # axs[1].set_ylim(0,1.25)
    axs[1,0].set_title('Normalized CF EPSCs at 10 Hz')
    # axs[1].set_xlabel('Pulse #')
    axs[1,0].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[1].legend(loc = "lower left")

    axs[2,0].errorbar(x, CF_WT_20hz_EPSCs, yerr = CF_WT_20hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[2,0].errorbar(x, CF_KO_20hz_EPSCs, yerr = CF_KO_20hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[2,0].errorbar(x, CF_WTCaBlock_20hz_EPSCs, yerr = CF_WTCaBlock_20hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[2,0].errorbar(x, CF_KOCaBlock_20hz_EPSCs, yerr = CF_KOCaBlock_20hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    # axs[2].set_xlim(1, 21)
    # axs[2].set_ylim(0,1.25)
    axs[2,0].set_title('Normalized CF EPSCs at 20 Hz')
    # axs[2].set_xlabel('Pulse #')
    axs[2,0].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[2].legend(loc = "lower left")

    axs[3,0].errorbar(x, CF_WT_50hz_EPSCs, yerr = CF_WT_50hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[3,0].errorbar(x, CF_KO_50hz_EPSCs, yerr = CF_KO_50hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[3,0].errorbar(x, CF_WTCaBlock_50hz_EPSCs, yerr = CF_WTCaBlock_50hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[3,0].errorbar(x, CF_KOCaBlock_50hz_EPSCs, yerr = CF_KOCaBlock_50hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    # axs[3].set_xlim(1, 21)
    # axs[3].set_ylim(0,1.25)
    axs[3,0].set_title('Normalized CF EPSCs at 50 Hz')
    axs[3,0].set_xlabel('Pulse #')
    axs[3,0].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[3].legend(loc = "lower left")

    #PLOT PF
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

    axs[0,1].errorbar(x, PF_WT_1hz_EPSCs, yerr = PF_WT_1hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[0,1].errorbar(x, PF_KO_1hz_EPSCs, yerr = PF_KO_1hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[0,1].errorbar(x, PF_WTCaBlock_1hz_EPSCs, yerr = PF_WTCaBlock_1hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[0,1].errorbar(x, PF_KOCaBlock_1hz_EPSCs, yerr = PF_KOCaBlock_1hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    #axs[0,1].set_xlim(1, 21)
    axs[0,1].set_ylim(0,4)
    axs[0,1].set_title('Normalized PF EPSCs at 1 Hz')
    #axs[0].set_xlabel('Pulse #')
    #axs[0,1].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[0,1].legend(loc = "lower left", ncol = 2)

    axs[1,1].errorbar(x, PF_WT_10hz_EPSCs, yerr = PF_WT_10hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[1,1].errorbar(x, PF_KO_10hz_EPSCs, yerr = PF_KO_10hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[1,1].errorbar(x, PF_WTCaBlock_10hz_EPSCs, yerr = PF_WTCaBlock_10hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[1,1].errorbar(x, PF_KOCaBlock_10hz_EPSCs, yerr = PF_KOCaBlock_10hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    # axs[1].set_xlim(1, 21)
    # axs[1].set_ylim(0,1.25)
    axs[1,1].set_title('Normalized PF EPSCs at 10 Hz')
    # axs[1].set_xlabel('Pulse #')
    #axs[1,1].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[1].legend(loc = "lower left")

    axs[2,1].errorbar(x, PF_WT_20hz_EPSCs, yerr = PF_WT_20hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[2,1].errorbar(x, PF_KO_20hz_EPSCs, yerr = PF_KO_20hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[2,1].errorbar(x, PF_WTCaBlock_20hz_EPSCs, yerr = PF_WTCaBlock_20hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[2,1].errorbar(x, PF_KOCaBlock_20hz_EPSCs, yerr = PF_KOCaBlock_20hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    # axs[2].set_xlim(1, 21)
    # axs[2].set_ylim(0,1.25)
    axs[2,1].set_title('Normalized PF EPSCs at 20 Hz')
    # axs[2].set_xlabel('Pulse #')
    #axs[2,1].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[2].legend(loc = "lower left")

    axs[3,1].errorbar(x, PF_WT_50hz_EPSCs, yerr = PF_WT_50hz_SEMs, marker = '.', linestyle='none', label = 'WT')
    axs[3,1].errorbar(x, PF_KO_50hz_EPSCs, yerr = PF_KO_50hz_SEMs, marker = '.', linestyle='none', label = 'KO')
    axs[3,1].errorbar(x, PF_WTCaBlock_50hz_EPSCs, yerr = PF_WTCaBlock_50hz_SEMs, marker = '.', linestyle='none', label = 'WT + Ca block')
    axs[3,1].errorbar(x, PF_KOCaBlock_50hz_EPSCs, yerr = PF_KOCaBlock_50hz_SEMs, marker = '.', linestyle='none', label = 'KO + Ca block')

    # axs[3].set_xlim(1, 21)
    # axs[3].set_ylim(0,1.25)
    axs[3,1].set_title('Normalized PF EPSCs at 50 Hz')
    axs[3,1].set_xlabel('Pulse #')
    #axs[3,1].set_ylabel(r'$EPSC_{n}/EPSC_{1}$')
    #axs[3].legend(loc = "lower left")

    plt.show()
