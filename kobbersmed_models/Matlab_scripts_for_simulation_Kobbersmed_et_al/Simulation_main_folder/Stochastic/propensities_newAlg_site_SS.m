function A_matrix_site = propensities_newAlg_site_SS(v_states, par, Ca_R_temp, SS_PM, Ca_prim_type) 

%Reaction rates of "site" (priming or calcium binding to second sensor)

m_max = par(14);
    
n = mod(v_states, 10);
m = mod(v_states-n,100)/10;

%Fast and slow sensor
k_4 = par(6);
k_min4 = par(7);
b_s = par(8); 
u_ves = par(43);

k_rep = par(10);
prim_kM = par(31);
prim_rate_const = par(32);
unprim_rate_const = par(33);
unprim_rate_const_0 = par(44);
    
if Ca_prim_type > 0
    [Caprim_rate, ~] = calculate_Caprim_rate(Ca_R_temp, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);
elseif u_ves~= 1
    Caprim_rate = prim_rate_const;
else
    Caprim_rate = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_matrix_site(1:3) = [  
        SS_PM*m*k_min4*(b_s).^((b_s~=0)*(m-1))...     %SmFn->S(m-1)Fn
        k_rep + Caprim_rate...        %Replenishment rate
        SS_PM*(m_max-m)*Ca_R_temp*k_4]; %SmFn->S(m+1)Fn

    
end