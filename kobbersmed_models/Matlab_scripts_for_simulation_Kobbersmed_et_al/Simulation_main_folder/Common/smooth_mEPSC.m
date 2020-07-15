function [mEPSC_time, mEPSC_current] = smooth_mEPSC

%Provides the mEPSC used in simulations

A = -7.209251536449789e-06; 
B = 2.709256850482493e-09;  
t_0 = 0; %0.002961444215447;
tau_rf = 10.692783377261414;
tau_df = 0.001500129264510;
tau_ds = 0.002823055510748;


mEPSC_gen_func = @(t)(t>=t_0).*(A*(1-exp(-(t-t_0)/tau_rf)).*(B*exp(-(t-t_0)/tau_df) + (1-B).*exp(-(t-t_0)/tau_ds)));

mEPSC_time = 0:1e-6:34*1e-3;

mEPSC_current = mEPSC_gen_func(mEPSC_time);