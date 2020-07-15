function test_run()

%%THIS SCRIPTS RUNS THE MODEL AT EXPERIMENTAL CALCIUM CONCENTRATIONS WITH
%%THE BEST FIT PARAMETERS OF THE UNPRIMING MODEL

model_type = 46; %The choice of model (and which parameters to define) as defined and explained in parameter_choices.m

%par_free are the choices of the parameters for hte choice of model_type
par_free(1) = 13.7718;
par_free(2) = 5;
par_free(3) = 5.5207e-08;
par_free(4) = 134.8542;
par_free(5) = 236.8238;
par_free(6) = 180;

stoch_on_off = 2; %0 for deterministic simulation, 1 for stochastic, 2 for both
rand_ves_on_off = 1; %This defines the vesicle placement as described in determ_vesicle_distances.m
CalC_on_off = 1; %0: No CalC simulation (if Ca files are already generated), >0 for CalC simulation. If == 1 all older calcium files are deleted
CaExtracellular = [0.75 1.5 3 6 10]; %Extracellular calcium concentrations (mM) to simulate
save_data = 2; %Defines how much data to be saved (1: A lot, including all eEJCs. 2: Data summary). See end of simulation_call_det.m and simulation_call_stoch.m
save_calc_loc = 1; %%Defines the location of calcium files (to be able to run more procedures in parallel without interfering with each other). If ==1, Calc simulation is 20.5 ms, if ==2 163.5 ms
pVr2_hack = 0; %Set this to zero. Only used for estimation of pVr2 (by putting new vesicles into the system right before second stimulation)

[par_init, savefilename] = parameter_choices(par_free, model_type, 0, rand_ves_on_off); %Defines the parameter vector and result filename to be inputted in the following function call. 

testing_the_system(stoch_on_off, rand_ves_on_off, CalC_on_off, par_init, CaExtracellular, save_data, savefilename, save_calc_loc, pVr2_hack)
