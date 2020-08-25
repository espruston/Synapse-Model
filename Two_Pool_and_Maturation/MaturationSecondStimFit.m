k_docking = .05;
k_undocking = 0;
k_maturation = .005;
k_dematuration = 0;

best = [0 0 0 0 1e30];
hz_compare = [.875; .577; .422; .319];
maturation = linspace(0.001,0.01,10);
ISIs = [1000 100 50 20];

for i = 1:10
    for j = 1:10
        for k = 1:10

            k_maturation = maturation(i);
            p_immature = .05*j;
            p_mature = .5+.05*k;

            delta_t = 1e-2; %ms

            t_SS = 10000; %ms
            ts_SS = linspace(0, t_SS, t_SS*delta_t + 1);

            state_0 = [0; 1; 0]; %start all vesicles in immature docked state 

            [t0,SS] = ode15s(@(t,state) dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration), [0 t_SS], state_0);

            ts = 0;
            Fused = zeros(2,1);
            second_stim = zeros(4,1);

            for m = 1:4
                stimulus_times = [0, ISIs(m)];

                max_time = stimulus_times(end) + stimulus_times(2)*3;
                stim_delay = diff(stimulus_times);
                stim_delay = [stim_delay max_time-stim_delay(end)];
                state = SS(end,:);

                for n = 1:2
                    pre_stim = state(end,:);
                    post_stim = pre_stim + [pre_stim(2)*p_immature+pre_stim(3)*p_mature -pre_stim(2)*p_immature -pre_stim(3)*p_mature];
                    Fused(n) = pre_stim(2)*p_immature + pre_stim(3)*p_mature;
                    [t,out] = ode45(@(t,state) dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration), [0 stim_delay(n)], post_stim);

                    state = [state(1:end-1,:); out];

                    ts = [ts(1:end-1); t+ts(end)];
                end
                second_stim(m) = Fused(2)/Fused(1);
            end

    %                 ts = [-10; -delta_t; ts];
    %                 state = [SS(end,:); SS(end,:); state];
            test_val = abs(second_stim - hz_compare);
            if mean(test_val)<best(end)
                best = [k_docking k_maturation p_immature p_mature mean(test_val)];
            end
        end
    end
end

function dydt = dSS(t,state,k_docking,k_undocking,k_maturation,k_dematuration)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end

function dydt = dState(t,state,k_docking,k_undocking,k_maturation,k_dematuration)
    
    dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
    dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
    dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;

end