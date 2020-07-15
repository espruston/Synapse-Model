function [CaRest] = calculate_CaRest(k_M_rest, CaMax, CaExt)

%%This script calculates resting calcium in M from given k_M and CaMax values
%Inputs in M.


CaRest = CaMax*(CaExt./(CaExt + k_M_rest));