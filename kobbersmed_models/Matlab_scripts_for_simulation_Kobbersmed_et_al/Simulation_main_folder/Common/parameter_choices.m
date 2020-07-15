function [par_init, savefilename] = parameter_choices(par_free, model_type, iteration, rand_ves_on_off)

%This file defines parameters of the simulation according to the choice of 'model_type'.
%The script also tests that the number of parameters match the model
%choice. Also provides an initial result file name (modified in other scripts).

%%Variables
%'par_free': Specific parameters for the chosen 'model_type'
%model_type: Choice of model. See list in code below.
%'iteration' is used when calling many repetitions (sets of 100).
%'rand_ves_on_off': only used when ==1000. It sets the
%number of vesicles to 1e12 for a good estimate of pVr2.


q = 1;

%%Model 1x: Single-sensor models.
%Model 11: par_free = Q_max
model_type_vec(q) = 11; num_pars_vec(q) = 1; q = q+1;
%Model 12: par_free = [Q_max kM]
model_type_vec(q) = 12; num_pars_vec(q) = 2; q = q+1;
%Model 13: par_free = [Q_max], k_rep = 0
model_type_vec(q) = 13; num_pars_vec(q) = 1; q = q+1;
%Model 14: par_free = [Q_max, num_ves], k_rep = 181.82
model_type_vec(q) = 14; num_pars_vec(q) = 2; q = q+1;
%Model 15: par_free = [Q_max, k_rep]
model_type_vec(q) = 15; num_pars_vec(q) = 2; q = q+1;
%Model 16: par_free = [Q_max, k_rep, num_ves]
model_type_vec(q) = 16; num_pars_vec(q) = 3; q = q+1;

%%Model 2x: Dual-sensor models. 
%Model 21: par_free = [Q_max, SS_coop, s]. SS_PM = 0.
model_type_vec(q) = 21; num_pars_vec(q) = 3; q = q+1;
%Model 22: par_free = [Q_max, SS_coop, s]. SS_PM = 1.
model_type_vec(q) = 22; num_pars_vec(q) = 3; q = q+1;
%Model 23: par_free = [Q_max, SS_coop, k_d_second, k_4, s]. SS_PM = 0
model_type_vec(q) = 23; num_pars_vec(q) = 5; q = q+1;
%Model 24: par_free = [Q_max, SS_coop, k_d_second, k_4, s]. SS_PM = 1
model_type_vec(q) = 24; num_pars_vec(q) = 5; q = q+1;
%Model 231: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep]. SS_PM = 0
model_type_vec(q) = 231; num_pars_vec(q) = 6; q = q+1;
%Model 232: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep, num_ves]. SS_PM = 0
model_type_vec(q) = 232; num_pars_vec(q) = 7; q = q+1;
%Model 241: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep]. SS_PM = 1
model_type_vec(q) = 241; num_pars_vec(q) = 6; q = q+1;
%Model 242: par_free = [Q_max, SS_coop, k_d_second, k_4, s, k_rep, num_ves]. SS_PM = 1
model_type_vec(q) = 242; num_pars_vec(q) = 7; q = q+1;
%Model 25: par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 0
model_type_vec(q) = 25; num_pars_vec(q) = 4; q = q+1;
%Model 26: par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 1
model_type_vec(q) = 26; num_pars_vec(q) = 4; q = q+1;
%Model 27: par_free = [Q_max, SS_coop, s]. SS_PM = 0, k_rep = 0;
model_type_vec(q) = 27; num_pars_vec(q) = 3; q = q+1;
%Model 28: par_free = [Q_max, SS_coop, s]. SS_PM = 1, k_rep = 0;
model_type_vec(q) = 28; num_pars_vec(q) = 3; q = q+1;

%%Model 3x: Site activation models
%Model 31: par_free = [Q_max, act_model_type, beta, gamma, delta]
model_type_vec(q) = 31; num_pars_vec(q) = 5; q = q+1;
%Model 32: par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep]
model_type_vec(q) = 32; num_pars_vec(q) = 6; q = q+1;
%Model 33: par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep, num_ves]
model_type_vec(q) = 33; num_pars_vec(q) = 7; q = q+1;

%%Model 4x: Unpriming models
%Model 41: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
model_type_vec(q) = 41; num_pars_vec(q) = 5; q = q+1;
%Model 42: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
model_type_vec(q) = 42; num_pars_vec(q) = 6; q = q+1;
%Model 43: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, k_M_Rest, CaMax_rest]
model_type_vec(q) = 43; num_pars_vec(q) = 7; q = q+1;
%Model 44: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0]
model_type_vec(q) = 44; num_pars_vec(q) = 6; q = q+1;
%Model 45: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, CaMax_rest]
model_type_vec(q) = 45; num_pars_vec(q) = 7; q = q+1;
%Model 46: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
model_type_vec(q) = 46; num_pars_vec(q) = 6; q = q+1;
%Model 47: par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
model_type_vec(q) = 47; num_pars_vec(q) = 5; q = q+1;

%Model 5x: Alternative unpriming model. NOT USED.
%Model 51: par_free = [Q_max, SS_coop, prim_rate, unprim_rate, u]; 
model_type_vec(q) = 51; num_pars_vec(q) = 5; %q = q+1;


%%Test that number of parameters match model type
model_ind = find(model_type_vec == model_type);
num_pars = num_pars_vec(model_ind);

if (length(par_free) ~= num_pars)
    error('Number of parameters is not correct for the model type')
end


%Basic assumptions: These may be changed if value in par_init specifies.

k_rep = 181.8182; %k_rep value if not specified in par_init

AZ_size = 0.6239936794; %Radius of simulation cylinder volume
kM = 2.679; %k_M of calcium influx ad resting calcium
CalC_geometry = 1; %1: Cylinder. 2: Box 
k_M_rest = kM;
CaMax_rest = 190e-9; %Aymptotic value of resting calcium (with increasing CaExt)

%Model 2
SS_coop = 0; %Coop. of second sensor 
SS_PM = 0; %Second sensor in plasma membrane or not
k_d_second = 0; %k_D of second sensor
k_on_second = 0; %k_on of second sensor
s_val = 1; %Fusion factor of second sensor
num_ves = 180; %Number of vesicles/sites
%Model 3
act_model_type = 0; %Cooperativity of site activation
act_const = 0; %Rate constant of activation
k_M_act = 0; %k_M of activation
gamma = 0; %Rate from [D] to [A]
delta = 0; %Rate from [A] to [D]
%Model 4
Ca_prim_type = 0; %Coop. of unpriming
prim_kM = 0; %k_M of unpriming
prim_rate_const = 0; %Priming rate
unprim_rate_const = 0; %Unpriming rate
unprim_rate_const_0 = 0; %Basal unpriming rate. Only if model_type == [44,45] (not used)
%Model 5
u_val = 1; %Factor on unpriming rate upon calcium binding (not used)

unprim_onestate = sum(model_type == [46 47]); %Only unpriming from state R[0].

if rand_ves_on_off == 1000
    num_ves = 1e12;
end


if sum(par_free ~= 0)
    %model_type 11: Model 1 with Q_max as an input. 
    %par_free = Q_max
    if model_type == 11
        Q_max = par_free(1);

    %model_type 12: Model 1 with Q_max and k_M as an input. 
    %par_free = [Q_max, kM]
    elseif model_type == 12
        Q_max = par_free(1);
        kM = par_free(2);

    %model_type 13: Model 1 with Q_max as an input. NO REPLENISHMENT. 
    %par_free = [Q_max]
    elseif model_type == 13
        Q_max = par_free(1);
        k_rep = 0;

     %model_type 14: Model 1 with Q_max and num_ves as an input. 
     %par_free = [Q_max, num_ves]
    elseif model_type == 14
        Q_max = par_free(1);
        num_ves = par_free(2);
        
    %model_type 15: Model 1 with Q_max and k_rep as an input. 
     %par_free = [Q_max, k_rep]
    elseif model_type == 15
        Q_max = par_free(1);
        k_rep = par_free(2);

    %model_type 15: Model 1 with Q_max and k_rep as an input. 
     %par_free = [Q_max, k_rep]
    elseif model_type == 16
        Q_max = par_free(1);
        k_rep = par_free(2);
        num_ves = par_free(3);

    %model_type 21: Model 2. 
    %par_free = [Q_max, SS_coop, s]. SS_PM = 0
    elseif model_type == 21
        Q_max = par_free(1);
        SS_coop = par_free(2);
        s_val = par_free(3);
        SS_PM = 0;
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
        
    %model_type 22: Model 2. 
    %par_free = [Q_max, SS_coop, s]. SS_PM = 1
    elseif model_type == 22
        Q_max = par_free(1);
        SS_coop = par_free(2);
        s_val = par_free(3);
        SS_PM = 1;
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
                
        %model_type 23: Model 2. 
        %par_free = [Q_max, SS_coop, k_4, k_min4, s]. SS_PM = 0
    elseif model_type == 23
        SS_PM = 0;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        k_d_second = par_free(3);
        k_on_second = par_free(4);
        s_val = par_free(5);
        
    %model_type 24: Model 2. 
    %par_free = [Q_max, SS_coop, k_4, k_min4, s]. SS_PM = 1
    elseif model_type == 24
        SS_PM = 1;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        k_d_second = par_free(3);
        k_on_second = par_free(4);
        s_val = par_free(5);

        %model_type 231: Model 2. 
        %par_free = [Q_max, SS_coop, k_4, k_min4, s, k_rep]. SS_PM = 0
    elseif model_type == 231
        SS_PM = 0;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        k_d_second = par_free(3);
        k_on_second = par_free(4);
        s_val = par_free(5);
        k_rep = par_free(6);

        %model_type 232: Model 2. 
        %par_free = [Q_max, SS_coop, k_4, k_min4, s, k_rep]. SS_PM = 0
    elseif model_type == 232
        SS_PM = 0;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        k_d_second = par_free(3);
        k_on_second = par_free(4);
        s_val = par_free(5);
        k_rep = par_free(6);
        num_ves = par_free(7);
        
    %model_type 241: Model 2. 
    %par_free = [Q_max, SS_coop, k_4, k_min4, s, k_rep]. SS_PM = 1
    elseif model_type == 241
        SS_PM = 1;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        k_d_second = par_free(3);
        k_on_second = par_free(4);
        s_val = par_free(5);
        k_rep = par_free(6);
    
    %model_type 242: Model 2. 
    %par_free = [Q_max, SS_coop, k_4, k_min4, s, k_rep, num_ves]. SS_PM = 1
    elseif model_type == 242
        SS_PM = 1;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        k_d_second = par_free(3);
        k_on_second = par_free(4);
        s_val = par_free(5);
        k_rep = par_free(6);
        num_ves = par_free(7);
        

    %model_type 25: Model 2. 
    %par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 0
    elseif model_type == 25
        Q_max = par_free(1);
        SS_coop = par_free(2);
        s_val = par_free(3);
        SS_PM = 0;
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
        num_ves = par_free(4);

    %model_type 26: Model 2. 
    %par_free = [Q_max, SS_coop, s, num_ves]. SS_PM = 1
    elseif model_type == 26
        Q_max = par_free(1);
        SS_coop = par_free(2);
        s_val = par_free(3);
        SS_PM = 1;
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
        num_ves = par_free(4);
        
    %model_type 27: Model 2. 
    %par_free = [Q_max, SS_coop, s]. SS_PM = 0
    elseif model_type == 27
        k_rep = 0;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        s_val = par_free(3);
        SS_PM = 0;
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
        
    %model_type 28: Model 2. 
    %par_free = [Q_max, SS_coop, s]. SS_PM = 1
    elseif model_type == 28
        k_rep = 0;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        s_val = par_free(3);
        SS_PM = 1;
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
                
        
    %model_type 31 Model 3. 
    %par_free = [Q_max, act_model_type, beta, gamma, delta]
    elseif model_type == 31
        act_const = 1e6;
        Q_max = par_free(1);
        act_model_type = par_free(2);
        beta = par_free(3);
        k_M_act = beta/act_const;
        gamma = par_free(4);
        delta = par_free(5);
        k_rep = par_free(6);
        
    %model_type 32 Model 3. 
    %par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep]
    elseif model_type == 32
        act_const = 1e6;
        Q_max = par_free(1);
        act_model_type = par_free(2);
        beta = par_free(3);
        k_M_act = beta/act_const;
        gamma = par_free(4);
        delta = par_free(5);
        k_rep = par_free(6);
        
    %model_type 33 Model 3. 
    %par_free = [Q_max, act_model_type, beta, gamma, delta, k_rep, num_ves]
    elseif model_type == 33
        act_const = 1e6;
        Q_max = par_free(1);
        act_model_type = par_free(2);
        beta = par_free(3);
        k_M_act = beta/act_const;
        gamma = par_free(4);
        delta = par_free(5);
        k_rep = par_free(6);
        num_ves = par_free(7);
        
    %model_type 41: Model 4. 
    %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
    elseif model_type == 41
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        
        %model_type 42: Model 4. 
        %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
    elseif model_type == 42
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        num_ves = par_free(6);
        
        %model_type 43: Model 4. 
        %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, k_M_Rest, CaMax_rest]
    elseif model_type == 43
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        k_M_rest = par_free(6);
        CaMax_rest = par_free(7);
        
        %Model_type 44: Model 4 with a basal unprim rate. 
        %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0]
    elseif model_type == 44
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        unprim_rate_const_0 = par_free(6);

        %Model_type 45: Model 4 with a basal unprim rate. 
        %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0]
    elseif model_type == 45
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        unprim_rate_const_0 = par_free(6);
        CaMax_rest = par_free(7);
        
        %model_type 46: Model 4. ONLY UNPRIMING FROM STATE [0,0]
        %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
    elseif model_type == 46
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        num_ves = par_free(6);

        %model_type 47: Model 4. ONLY UNPRIMING FROM STATE [0,0]
        %par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
    elseif model_type == 47
        k_rep = 0;
        Q_max = par_free(1);
        Ca_prim_type = par_free(2);
        prim_kM = par_free(3);
        prim_rate_const = par_free(4);
        unprim_rate_const = par_free(5);
        
    elseif model_type == 51
        k_rep = 0;
        Q_max = par_free(1);
        SS_coop = par_free(2);
        prim_rate_const = par_free(3);
        unprim_rate_const = par_free(4);
        u_val = par_free(5);
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;

    end
    
elseif par_free == 0    
    %model_type 11: Model 1 with Q_max as an input. par_free = Q_max
    if model_type == 11
        [par_free(1), Q_max] = deal(5); %5.693;

    %model_type 12: Model 1 with Q_max and k_M as an input. par_free = [Q_max, kM]
    elseif model_type == 12
        [par_free(1), Q_max] = deal(5); %5.693;
        [par_free(2), kM] = deal(2.679);

    %model_type 13: Model 1 with Q_max as an input. NO REPLENISHMENT. par_free = [Q_max]
    elseif model_type == 13
        [par_free(1), Q_max] = deal(5); %5.693;
        [par_free(2), k_rep] = deal(0);

    %model_type 21/23: Model 2. par_free = [Q_max, SS_coop, k_4, k_min4,
    %s]. SS_PM = 0.
    elseif model_type == 21 || model_type == 23
        [par_free(1), Q_max] = deal(5);
        [par_free(2), SS_coop] = deal(2);
        [par_free(3), SS_PM] = deal(0);
        [par_free(4), k_on_second] = deal((1.4e8)/2);
        [par_free(5), k_off_second] = deal(4000/60);
        [par_free(6), k_d_second] = deal(k_off_second/k_on_second);
        [par_free(7), s_val] = deal(45);
        
        %model_type 22/24: Model 2. par_free = [Q_max, SS_coop, k_4,
        %k_min4, s]. SS_PM = 1
    elseif model_type == 22 || model_type == 24
        [par_free(1), Q_max] = deal(5);
        [par_free(2), SS_coop] = deal(2);
        [par_free(3), SS_PM] = deal(1);
        [par_free(4), k_on_second] = deal((1.4e8)/2);
        [par_free(5), k_off_second] = deal(4000/60);
        [par_free(6), k_d_second] = deal(k_off_second/k_on_second);
        [par_free(7), s_val] = deal(45);
    
    %model_type 31 Model 3. par_free = [Q_max, act_model_type, kMact, gamma, delta]
    elseif model_type == 31
        [par_free(1), Q_max] = deal(12);
        [par_free(2), act_model_type] = deal(3);
        [par_free(3), k_M_act] = deal(120e-9);
        [par_free(4), gamma] = deal(150);
        [par_free(5), delta] = deal(10);

    %model_type 41: Model 4. par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const]
    elseif model_type == 41
        [k_rep] = deal(0);
        [par_free(1), Q_max] = deal(12);
        [par_free(2), Ca_prim_type] = deal(3);
        [par_free(3), prim_kM] = deal(55e-9);
        [par_free(4), prim_rate_const] = deal(320);
        [par_free(5), unprim_rate_const] = deal(200);
        
    %model_type 42: Model 4. par_free = [Q_max, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, num_ves]
    elseif model_type == 42
        [k_rep] = deal(0);
        [par_free(1), Q_max] = deal(12);
        [par_free(2), Ca_prim_type] = deal(3);
        [par_free(3), prim_kM] = deal(55e-9);
        [par_free(4), prim_rate_const] = deal(320);
        [par_free(5), unprim_rate_const] = deal(200);
        [par_free(6), num_ves] = deal(180);
        
    elseif model_type == 43
        [k_rep] = deal(0);
        [par_free(1), Q_max] = deal(12);
        [par_free(2), Ca_prim_type] = deal(3);
        [par_free(3), prim_kM] = deal(55e-9);
        [par_free(4), prim_rate_const] = deal(320);
        [par_free(5), unprim_rate_const] = deal(200);
        [par_free(6), k_M_rest] = deal(2.697);
        [par_free(7), CaMax_rest] = deal(190e-9);
        
    elseif model_type == 44
        [k_rep] = deal(0);
        [par_free(1), Q_max] = deal(12);
        [par_free(2), Ca_prim_type] = deal(3);
        [par_free(3), prim_kM] = deal(55e-9);
        [par_free(4), prim_rate_const] = deal(320);
        [par_free(5), unprim_rate_const] = deal(200);
        [par_free(6), unprim_rate_const_0] = deal(30);
        
    elseif model_type == 51
        k_rep = 0;
        [par_free(1), Q_max] = deal(12);
        [par_free(2), SS_coop] = deal(2);
        [par_free(3), prim_rate_const] = deal(320);
        [par_free(4), unprim_rate_const] = deal(200);
        [par_free(5), u_val] = deal(100);
        k_on_second = (1.4e8)/2;
        k_off_second = 4000/60;
        k_d_second = k_off_second/k_on_second;
    end    
end





%% par_init is defined below

%Basic
par_init(1) = CalC_geometry;
par_init(2) = Q_max;
par_init(3) = kM;
par_init(4) = num_ves;
par_init(5) = k_rep;
%Model 2
par_init(6) = SS_coop;
par_init(7) = SS_PM;
par_init(8) = k_d_second;
par_init(9) = k_on_second;
par_init(10) = s_val;
%Model 3
par_init(11) = act_model_type;
par_init(12) = k_M_act;
par_init(13) = gamma;
par_init(14) = delta;
%Model 4
par_init(15) = Ca_prim_type;
par_init(16) = prim_kM;
par_init(17) = prim_rate_const;
par_init(18) = unprim_rate_const;

par_init(19) = k_M_rest;
par_init(20) = CaMax_rest;

par_init(21) = unprim_rate_const_0;

%Model 5
par_init(22) = u_val;


%Model 3, extra:
par_init(23) = act_const;

%Common, extra
par_init(24) = AZ_size;

par_init(25) = model_type;

par_init(26) = unprim_onestate;



filesuf = 'parfree';

for l = 1:length(par_free)
    filesuf = [filesuf '_' num2str(par_free(l))];
end

CaMax_filename = CaMax_rest*1e9;

if iteration == 0 %In deterministic simulations and stochastic with 200 repetitions.
    savefilename = ['Results_CaMax' num2str(CaMax_filename) '_model' num2str(model_type) filesuf];
else %If running many repetitions ('run_more_reps') several result files are generated and concenated. 
    savefilename = ['Results' num2str(iteration) '_CaMax' num2str(CaMax_filename) '_model' num2str(model_type) filesuf];
end


