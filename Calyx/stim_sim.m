function [ts, state, Fused_im, Fused_m, Ca_sim] = stim_sim(stimulus_times, max_time)

    global p_immature
    global p_mature
    global k_docking
    global k_undocking
    global k_maturation
    global k_dematuration
    global k_on_3
    global k_off_3
    global k_on_7
    global k_off_7
    global C_3
    global C_7
    global delta_t
    global SS
    
    state = SS;
    Ca_sim = create_Ca_signal(stimulus_times, max_time);
    %[syt3, syt7] = syt_sim(Ca_sim);
    stim_delay = diff(stimulus_times);
    stim_delay = [stim_delay max_time-stimulus_times(end)];

    ts = 0;
    Fused_im = zeros(length(stimulus_times),1);
    Fused_m = zeros(length(stimulus_times),1);

    for i = 1:length(stim_delay)

        pre_stim = state(end,:);
        post_stim = pre_stim + [pre_stim(2)*p_immature+pre_stim(3)*p_mature, -pre_stim(2)*p_immature, -pre_stim(3)*p_mature, 0, 0];
        Fused_im(i) = pre_stim(2)*p_immature;
        Fused_m(i) = pre_stim(3)*p_mature;
        [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration,k_on_3,k_off_3,k_on_7,k_off_7,C_3,Ca_sim), [ts(end) ts(end)+stim_delay(i)], post_stim);

        state = [state(1:end-1,:); out];

        ts = [ts(1:end-1,:); t];
    end

    state = [SS; state];
    ts = [delta_t; ts];
end