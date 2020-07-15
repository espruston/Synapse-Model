function call_result_check_new

CalC_on_off = 1;

save_data = 1;
stoch_on_off = 2;
rand_ves_on_off = 1;
par_free = 0;

model_type = 11

[par_init] = parameter_choices(par_free, model_type);

testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data);

model_type = 22
CalC_on_off = 0;
[par_init] = parameter_choices(par_free, model_type);

testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data);

model_type = 23

[par_init] = parameter_choices(par_free, model_type);

testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data);

tic;

model_type = 31
CalC_on_off = 1;
[par_init] = parameter_choices(par_free, model_type);

testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data);
toc

tic;
model_type = 41
CalC_on_off = 0;
[par_init] = parameter_choices(par_free, model_type);

testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, save_data);
toc