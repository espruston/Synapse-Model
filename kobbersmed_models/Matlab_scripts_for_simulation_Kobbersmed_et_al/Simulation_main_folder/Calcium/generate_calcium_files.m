function generate_calcium_files()


rand_ves_on_off = 0;
Q_maxes = 4.5169;

parfor k = 1:length(Q_maxes)*5
    Q_ind = ceil(k/5);
    par_free = Q_maxes(Q_ind);
    
    cal_ind = mod(k-1,5)+1;
    
    [par_init, ~] = parameter_choices(par_free, 11, 0, 0)

   
    par = model_parameters_det(par_init);
    CaExtracellular = [0.75 1.5 3 6 10];


    % %%%%%Calcium, residual and vesicles%%%%
    CaMax_rest = par(23); 
    k_M_rest = par(24);
    CaRest = calculate_CaRest(k_M_rest, CaMax_rest, CaExtracellular);

    long_sim_on_off = 1;


    RunCalC_det_altsaving(CaExtracellular(cal_ind), par, CaRest(cal_ind), long_sim_on_off); %Bm_conc, ATP_conc, fast_buffer_on_off, fast_buffer_params, 

end