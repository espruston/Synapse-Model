function Ca_sim = create_Ca_signal(stimulus_times, max_time)
    
    global Ca_rest
    global Ca_spike
    global Ca_residual
    global T_Ca_decay
    global delta_t
    global sigma
    global mu
    
    
    ts = linspace(0,max_time,max_time/delta_t + 1);
    Ca_sim = zeros(1,length(ts));
    Ca_sim = Ca_sim + Ca_rest;

    for t = 1:length(stimulus_times) %simulate calcium influx

        spike_start_index = round(stimulus_times(t)/delta_t) + 1; %if 1st stim is at t=0 index should be one, round is necessary due to IEEE fp returning scientific notation ocasionally
        spike_peak_index = round((stimulus_times(t)+mu)/delta_t) + 1;

        Ca_sim(spike_start_index:end) = Ca_sim(spike_start_index:end) + Ca_spike*exp(-1*((ts(1:end - spike_start_index + 1) - mu)/sigma).^2/2); %if 1st stim is at t=0 index should be one

        Ca_sim(spike_peak_index:end) = Ca_sim(spike_peak_index:end) + Ca_residual*exp(-1*ts(1:end - spike_peak_index + 1)/T_Ca_decay);    

    end

end