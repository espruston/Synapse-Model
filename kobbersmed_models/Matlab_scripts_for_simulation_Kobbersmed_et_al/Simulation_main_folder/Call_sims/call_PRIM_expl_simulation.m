function call_PRIM_expl_simulation()

Ca_prim_type = 3;%0;
% prim_kMs = [60e-9 80e-9 100e-9 120e-9];%0;
prim_kMs = [40e-9 55e-9 65e-9 70e-9 75e-9 100e-9 120e-9 150e-9 200e-9];
prim_rate_consts = [100 150 180 210 250 300 500];%0;
unprim_rate_consts = [50 100 150 200 300 400 500];%0;

Q_max = 16;

model_type = 41;
stoch_on_off = 0;
rand_ves_on_off = 1;
save_data = 2;

% for l = 1:2
%     SS_PM = SS_PM_vals(l)


%         if k > 1
%            CalC_on_off = 0; 
%         end

CalC_on_off = 1;
par_free = [Q_max Ca_prim_type 0 0 0];

[par_init, savefilename] = parameter_choices(par_free, model_type)

tic; testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
toc



CalC_on_off = 0;
for k = 1:length(prim_kMs)
    prim_kM = prim_kMs(k)
    for l = 1:length(prim_rate_consts)
        prim_rate_const = prim_rate_consts(l)
        for m = 1:length(unprim_rate_consts)
            

            
            
            unprim_rate_const = unprim_rate_consts(m)
            
            par_free = [Q_max Ca_prim_type prim_kM prim_rate_const unprim_rate_const];
                
            [par_init, savefilename] = parameter_choices(par_free, model_type)
            
            tic; testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data, savefilename);
            toc
        end
    end
end

