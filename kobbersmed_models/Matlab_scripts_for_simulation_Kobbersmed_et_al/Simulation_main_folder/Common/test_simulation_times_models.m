function test_simulation_times_models()

save_data = 1;
rand_ves_on_off = 1;
% par_free = 0;


model_types = [11 22 23 31 41];

model_time = struct('type', cell(1, length(model_types)*2), 'simtime', cell(1, length(model_types)*2));

for k = 1:(2*length(model_types))

    %Run CalC
    if k == 1
        CalC_on_off = 1;
        model_type = 11;
        par_free = 5;
        [par_init] = parameter_choices(par_free, model_type);
        testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, save_data);
    elseif k == 7
        CalC_on_off = 1;
        model_type = 11;
        par_free = 12;
        [par_init] = parameter_choices(par_free, model_type);
        testing_the_system(999, rand_ves_on_off, CalC_on_off, par_init, save_data);   
    end

    CalC_on_off = 0;
    
    model_ind = ceil(k/2);
    stoch_on_off = mod(k-1,2);
    model_type = model_types(model_ind);
    if stoch_on_off == 0
        model_time(k).type = ['model_type' num2str(model_type) 'DET'];
    elseif stoch_on_off == 1
        model_time(k).type = ['model_type' num2str(model_type) 'STOCH'];
    end
    
    [par_init] = parameter_choices(0, model_type);
    tic;
    testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data);
    model_time(k).simtime = toc
    
    save('./Sim_data/new_results/SimTimeTestResult', 'model_time')
end



