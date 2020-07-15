function [Caprim_rate, Caunprim_rate] = calculate_Caprim_rate(Calcium, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0)


Caprim_rate = (Ca_prim_type~=0)*prim_rate_const;

Caunprim_rate = (Ca_prim_type~=0)*(unprim_rate_const*(1 - (Calcium.^Ca_prim_type ./(Calcium.^Ca_prim_type + prim_kM))) + unprim_rate_const_0);