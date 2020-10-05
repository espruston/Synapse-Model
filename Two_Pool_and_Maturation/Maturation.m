classdef Maturation
    %Maturation is a class containing maturation-based models of synaptic
    %vesicle release
    %   Detailed explanation goes here
    
    properties
        
        p_immature
        p_mature
        state
        ts
        syts
        k_maturation
        k_dematuration
        k_docking
        k_undocking
        k_on_3
        k_off_3
        k_on_7
        k_off_7
        C_3
        C_7
        Ca_rest
        Ca_spike
        Ca_residual
        T_Ca_decay
        stimulus_times
        delta_t
        max_time
        Ca_sim
        t_SS
        SS
        
    end
    
    methods
        function obj = Maturation(p_immature,p_mature,k_maturation,k_dematuration,k_docking,k_undocking,k_on_3,k_off_3,k_on_7,k_off_7,C_3,C_7)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.p_immature = p_immature;
            obj.p_mature = p_mature;
            obj.k_maturation = k_maturation;
            obj.k_dematuration = k_dematuration;
            obj.k_docking = k_docking;
            obj.k_undocking = k_undocking;
            obj.k_on_3 = k_on_3;
            obj.k_off_3 = k_off_3;
            obj.k_on_7 = k_on_7;
            obj.k_off_7 = k_off_7;
            obj.C_3 = C_3;
            obj.C_7 = C_7;
            
            obj.t_SS = 10000;
            
        end
        
        function obj = stim_sim(obj,Ca_rest,Ca_spike,Ca_residual,T_Ca_decay,delta_t,sigma,mu,stimulus_times,max_time)
            
            obj.SS = steady_state(obj); 
            
            obj.Ca_rest = Ca_rest;
            obj.Ca_spike = Ca_spike;
            obj.Ca_residual = Ca_residual;
            obj.T_Ca_decay = T_Ca_decay;
            obj.delta_t = delta_t;
            obj.stimulus_times = stimulus_times;
            obj.max_time = max_time;
            obj.Ca_sim = create_Ca_signal(Ca_rest, Ca_spike, Ca_residual, T_Ca_decay, delta_t, sigma, mu, stimulus_times, max_time);
            
            state = obj.SS;
            stim_delay = diff(stimulus_times);
            stim_delay = [stim_delay max_time-stimulus_times(end)];

            ts = 0;
            Fused_im = zeros(length(stimulus_times),1);
            Fused_m = zeros(length(stimulus_times),1);

            for i = 1:length(stim_delay)

                pre_stim = state(end,:);
                post_stim = pre_stim + [pre_stim(2)*obj.p_immature*pre_stim(4)+pre_stim(3)*obj.p_mature, -pre_stim(2)*obj.p_immature*pre_stim(4), -pre_stim(3)*obj.p_mature];
                Fused_im(i) = pre_stim(2)*obj.p_immature*pre_stim(4);
                Fused_m(i) = pre_stim(3)*obj.p_mature;
                [t,out] = ode45(@(t,state) dState(t,state,obj.k_docking,obj.k_undocking,obj.k_maturation,obj.k_dematuration,obj.Ca_sim), [ts(end) ts(end)+stim_delay(i)], post_stim);

                state = [state(1:end-1,:); out];

                ts = [ts(1:end-1,:); t];
            end

            obj.state = [obj.SS; state];
            obj.ts = [delta_t; ts];
            
        end
        
        function SS = steady_state(obj)
            
            state_0 = [1;0;0]; %all sites are empty to start
            [t,State] = ode15s(@(t,State) dSS(t,State,obj.k_docking,obj.k_undocking,obj.k_maturation,obj.k_dematuration), [0 obj.t_SS], state_0);

            SS = State(end,:);
            
        end
        
        function obj = syt_sim(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            syts_0 = [0;0];
            [t, syts] = ode45(@(t,syts) dSyt(t,syts,obj.k_on_3,obj.k_off_3,obj.k_on_7,obj.k_off_7,obj.Ca_sim,obj.delta_t), [0, obj.max_time], syts_0);
            obj.syts = syts;
        end
        
        function dydt = dSyt(t,syts,k_on_3,k_off_3,k_on_7,k_off_7,Ca_sim,delta_t)
            Ca = Ca_sim(round(t/delta_t)+1);
            
            dydt(1,1) =  (1-syts(1))*k_on_3*Ca - syts(1)*k_off_3;
            dydt(2,1) =  (1-syts(2))*k_on_7*Ca - syts(2)*k_off_7;
            
        end
        
        function dydt = dSS(~,state,k_docking,k_undocking,k_maturation,k_dematuration)
            
            dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
            dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
            dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
            
        end
        
        function dydt = dState(~,state,k_docking,k_undocking,k_maturation,k_dematuration)
            
            dydt(1,1) = -state(1)*k_docking + state(2)*k_undocking;
            dydt(2,1) = state(1)*k_docking - state(2)*k_undocking - state(2)*k_maturation + state(3)*k_dematuration;
            dydt(3,1) = state(2)*k_maturation - state(3)*k_dematuration;
            
        end
    end
end

