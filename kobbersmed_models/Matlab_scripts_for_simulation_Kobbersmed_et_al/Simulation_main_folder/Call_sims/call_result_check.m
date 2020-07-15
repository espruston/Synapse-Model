function call_result_check()

large_simulation = 1;

Ca_prim_type = 3;%0;
prim_kM = 55e-9;%0; 60e-9
prim_rate_const = 320;%0; 180
unprim_rate_const = 200;%0; 300

k_rep = 0; %Ca-independent replenishment

num_ves = 180;% * num_ves_factor;
act_model_type = 0;
k_M_act_inp = 0; %120e-9;
gamma = 0;%0;
delta = 0;%0;
% Q_maxes = 0.5:0.5:12.5;

%%%One-sensor-model

k_3_first = 1.4e8;
kmin3 = 4000;
k_d_first = kmin3/k_3_first;

%Regehr values
SS_PM = 0;
SS_coop = 0;
opt_k_4_second = 0;%
kmin4 = 0;
opt_k_d_second = 0;
opt_s = 1;

Q_max = 12;

CalC_on_off = 1

rand_ves_on_off = 1;

% for l = 1:2
%     SS_PM = SS_PM_vals(l)


        if k > 1
           CalC_on_off = 0; 
        end


            tic; testing_the_system(2, act_model_type, Ca_prim_type, SS_PM, SS_coop, k_M_act_inp, rand_ves_on_off, CalC_on_off, Q_max, num_ves, k_rep, gamma, delta, prim_kM, prim_rate_const, unprim_rate_const, opt_s, opt_k_d_second, opt_k_4_second, 2);
 
            toc
% end








if large_simulation == 1

    CalC_on_off = 0
    num_ves = 180;
    height = 1;


    Ca_prim_type = 0;
    prim_kM = 0;
    prim_rate_const = 0;
    unprim_rate_const = 0;

    k_rep = 180; %Ca-independent replenishment

    SS_PM = 0;
    SS_coop = 0;
    opt_k_4_second = 0;%
    kmin4 = 0;
    opt_k_d_second = 0;
    opt_s = 1;

    Q_max = 12;
    act_model_type = 3;
    k_M_act_inp = 120e-9; %125e-9;
    gamma = 150;
    delta = 10;

%     for SS_PM = [0 1]
        tic; testing_the_system(2, act_model_type, Ca_prim_type, SS_PM, SS_coop,  k_M_act_inp, 1, CalC_on_off, Q_max, num_ves, k_rep, gamma, delta, prim_kM, prim_rate_const, unprim_rate_const, opt_s, opt_k_d_second, opt_k_4_second, 2);
        toc
%     end





    num_ves = 180;% * num_ves_factor;
    act_model_type = 0;
    k_M_act_inp = 0; %120e-9;
    gamma = 0;%0;
    delta = 0;%0;
    % Q_maxes = 0.5:0.5:12.5;
    height = 1;

    %%%One-sensor-model

    k_3_first = 1.4e8;
    kmin3 = 4000;
    k_d_first = kmin3/k_3_first;

    %Regehr values
    SS_PM_vals = [0 1];
    SS_coop = 2;
    opt_k_4_second = k_3_first/2;%  0;
    kmin4 = kmin3/60;
    opt_k_d_second = kmin4/opt_k_4_second;%   0;
    opt_s_vals = 45;

    Q_max = 3;

    CalC_on_off = 1

    rand_ves_on_off = 1;

    for l = 1:2
        SS_PM = SS_PM_vals(l)
            opt_s = opt_s_vals

            if l > 1
               CalC_on_off = 0; 
            end

        tic; testing_the_system(2, act_model_type, Ca_prim_type, SS_PM, SS_coop, k_M_act_inp, rand_ves_on_off, CalC_on_off, Q_max, num_ves, k_rep, gamma, delta, prim_kM, prim_rate_const, unprim_rate_const, opt_s, opt_k_d_second, opt_k_4_second, 2);
        toc
    end


    %%%One-sensor-model
    SS_PM = 0;
    SS_coop = 0;
    opt_k_4_second = 0;%  0;
    kmin4 = 0;
    opt_k_d_second = 0;%   0;
    opt_s = 1;

    Q_max = 3;

    CalC_on_off = 0

    rand_ves_on_off = 1;



        tic; testing_the_system(2, act_model_type, Ca_prim_type, SS_PM, SS_coop, k_M_act_inp, rand_ves_on_off, CalC_on_off, Q_max, num_ves, k_rep, gamma, delta, prim_kM, prim_rate_const, unprim_rate_const, opt_s, opt_k_d_second, opt_k_4_second, 2);
        toc

end







