function A_matrix = propensities_newAlg_SS(v_states, par, Ca_R_temp, fuse_state, act_number_vesicles, SS_PM, Ca_prim_type) %Returns the propensities of the PCCS-model reactions

%Reaction rates of the whole system to determine reaction time and reacting
%vesicle.
    
n_max = par(13);
m_max = par(14);
    
    
n = mod(v_states, 10);
m = mod((v_states-n),100)/10;

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
k_rep = par(10);
u_val = par(43);
unprim_onestate = par(47);

if f == 0 || s == 0
    l_0 = 0;
end


prim_kM = par(31);
prim_rate_const = par(32);
unprim_rate_const = par(33);
unprim_rate_const_0 = par(44);
[prim_rates, unprim_rates] = calculate_Caprim_rate(Ca_R_temp, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0);

if u_val~=1
    prim_rates = prim_rate_const;
    unprim_rates = (v_states<fuse_state).*unprim_rate_const.*u_val.^(-1*m);
end


if unprim_onestate
    unprim_rates = (v_states == 0).*unprim_rates;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



A_matrix = ...
    (v_states >= fuse_state) .* (SS_PM .* ((m_max-m).*Ca_R_temp*k_4 + m.*k_min4.*(b_s).^((b_s~=0)*(m-1))) + k_rep + prim_rates) ...
    + (v_states < fuse_state).*((act_number_vesicles == 1) .* l_0.*f.^n.*s.^m + (n_max-n).*Ca_R_temp*k_3 + n.*k_min3.*(b_f).^((b_f~=0)*(n-1)) + (m_max-m).*Ca_R_temp*k_4 + m.*k_min4.*(b_s).^((b_s~=0)*(m-1)) + unprim_rates) ...
                ;  


end