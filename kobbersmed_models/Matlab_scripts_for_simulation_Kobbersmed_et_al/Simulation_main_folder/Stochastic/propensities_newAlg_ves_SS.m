function [A_matrix_ves, n, m] = propensities_newAlg_ves_SS(ves_states, par, Ca_R_temp, act_status, Ca_prim_type) %Returns the propensities of the PCCS-model reactions

%Reaction rates of single SV

n_max = par(13);
m_max = par(14);
    
n = mod(ves_states, 10);
m = (ves_states-n)/10;

%Fast and slow sensor
k_3 = par(1);
k_min3 = par(2);
b_f = par(3);
f = par(5);
l_0 = par(4); 
k_4 = par(6);
k_min4 = par(7);
b_s = par(8);    
s = par(9);
u_val = par(43);
unprim_onestate = par(47);



if f == 0 || s == 0
    l_0 = 0;
end



unprim_rate_const = par(33);

if Ca_prim_type > 0
    prim_kM = par(31);
    prim_rate_const = par(32);
    unprim_rate_const_0 = par(44);
    [~, unprim_rate] = calculate_Caprim_rate(Ca_R_temp, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);

    if unprim_onestate
        unprim_rate = unprim_rate * (ves_states == 0);
    end
    
elseif u_val ~= 1
    unprim_rate = unprim_rate_const*u_val^(-1*m);
else
    unprim_rate = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_matrix_ves(1:6) = [  
        m*k_min4*(b_s).^((b_s~=0)*(m-1))...     %SmFn->S(m-1)Fn
    n*k_min3*(b_f).^((b_f~=0)*(n-1))...  %SmFn->SmF(n-1)    
    act_status * l_0*f^n.*s.^m...        %SmFn->Fuse
    (n_max-n)*Ca_R_temp*k_3 ...  %SmFn->SmF(n+1)
    (m_max-m)*Ca_R_temp*k_4 ... %SmFn->S(m+1)Fn
     unprim_rate];

    
end