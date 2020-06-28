import pandas as pd
import matplotlib as plt

if __name__ == "__main__":

    CF_data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skip_rows = 1, nrows = 20)
    CF_data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 23, nrows = 20)
    CF_data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 45, nrows = 20)
    CF_data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 4, usecols = 'B,C,F,G,J,K,N,O', skiprows = 67, nrows = 20)

    PF_data_1hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skip_rows = 1, nrows = 20)
    PF_data_10hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 23, nrows = 20)
    PF_data_20hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 45, nrows = 20)
    PF_data_50hz = pd.read_excel(r'~/Dropbox/Work/Jackman Lab/Modeling/Models and raw data_Dennis.xlsx', sheet_name = 6, usecols = 'B,C,F,G,J,K,N,O', skiprows = 67, nrows = 20)

    
