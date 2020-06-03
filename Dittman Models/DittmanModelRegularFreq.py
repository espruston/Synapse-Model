from math import e

#table 2 parameters
roh = 3.1
F_1 = .05
T_F = 100
T_D = 50
k_0 = 2/1000
k_max = 30/1000
K_D = 2

#other parameters
r = 50 #frequency in hz
T_E = 2
delta_F =
delta_D =
D_1 = 1

alpha =

#steady state values
CaX_F_ss = delta_F*(e**(1/(r*T_D)) - 1)**(-1)

CaX_D_ss = delta_D*(1-e**(-1/(r*T_D)))

F_ss = F_1 + (1-F_1)/(1+(K_F/CaX_F_ss))

xi_ss =

D_ss = (1 - e**(-1*k_0/r)*xi_ss)/(1 - (1 - F_ss)*e**(-1*k_0/r)*xi_ss)

EPSC_ss =

#indexed values
CaX_F[i] = CaX_F_ss*(1-e**(-1*(i-1)/(r*T_D)))

F[i] = F_1 + (1-F_1)/(1+(K_F/CaX_F[i]))

CaX_D[i-1] = CaX_D_ss*(1-e**(-1*(i-1)/(r*T_D)))

xi[i] = ((K_D/CaX_D[i-1] + 1)/((K_D/CaX_D[i-1]) + e**(-1/(r*T_D))))**(-1*(k_max-k_0)*T_D) #equation 16

D[i] = 1 - (1 - (1 - F[i-1]))*(D[i-1])*e**(k_0/r)*xi[i] #equation 15

EPSC[i] = N_T*D[i]*F[i]*(t*e/T_E)*e**(-t/T_E)

#equation 19
EPSC_norm[i] = EPSC[i]/EPSC[0]
EPSC_norm[i] = (N_T*D[i]*F[i])/(N_T*D_1*F_1)
EPSC_norm[i] = D[i]*(F[i]/F_1)

#equation 21
EPSC_norm_ss = EPSC_ss/EPSC[0]
EPSC_norm_ss = D_ss*(F_ss/F_1)
