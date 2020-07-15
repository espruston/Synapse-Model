function [ t, fused_ves, ves_states, act_states] = directMethod_newAlg_steps_SS(tspan, ves_init_states, ss_act_number, par,...
                                  Ca_time_vesicles, Ca_R_vesicles, n_max, m_max, fuse_state, SS_PM, num_ves, act_model_type, Ca_prim_type, collect_states) 
%DIRECTMETHOD Implementation of the Direct Method variant of the Gillespie algorithm
%   Based on: Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%   Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%%  Author: Janus R?nn Lind, 2016 <januslind@math.ku.dk>
%
%   This is a special case of the use of Gillespise algorithm

min_step = par(37);

if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1e6;
end
if ~exist('par', 'var')
    par = [];
end

if n_max > 9 || m_max > 9
    disp('Warning. Number of binding sites exceeds limit.')
    return
end
    
%% Initialize
T = zeros(MAX_OUTPUT_LENGTH, 1);
fused_ves = zeros(MAX_OUTPUT_LENGTH, 1);
T(1)     = tspan(1);
X   = ves_init_states;
rxn_count = 1; %Reaction counter
T_cur = 0;
tau_cur = 0;


%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%
if collect_states
    ves_states = zeros(length(ves_init_states), MAX_OUTPUT_LENGTH);
    act_states = zeros(length(ss_act_number), MAX_OUTPUT_LENGTH);

    ves_states(:, 1) = ves_init_states;
    act_states(:, 1) = ss_act_number;
end
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%

act_number_vesicles = ss_act_number';

act_rate_const = par(26);
inact_rate_const = par(27);

delay_rate = par(28);


%% MAIN LOOP
while T_cur + tau_cur < tspan(2)
    
%     if ((T(rxn_count) + tau_cur > 0.04) && (T(rxn_count) + tau_cur <0.01)) || (T(rxn_count) + tau_cur >0.014)
%         min_step = 1e-4;
%     else
%         min_step = min_step_init;
%     end
    
    % Determine calcium at relevant time point
    
    Ca_R_temp = interp1q(Ca_time_vesicles, Ca_R_vesicles, T_cur+tau_cur)'; %Current Ca concentration, vesicles
    
    % Calculate reaction propensities
    a1  = propensities_newAlg_SS(X, par, Ca_R_temp, fuse_state, act_number_vesicles, SS_PM, Ca_prim_type);
    

    a2 = propensities_newAlg_act(act_number_vesicles, par, Ca_R_temp, act_model_type);
    

    a = [a1; a2]; 
    
    % Compute tau and mu
    a0 = sum(a);
    
    if a0 > 0 %Condition makes sure there are still releaseable vesicles

        r = rand(1,3);  
        tau = (1/a0)*log(1/r(1));

        if tau > min_step
           tau_cur = tau_cur + min_step;
        else
            ves_ind  = find((cumsum(a) >= r(2)*a0),1,'first'); % Determine vesicle number

            if rxn_count + 1 > MAX_OUTPUT_LENGTH
                t = T(1:rxn_count);
                fused_ves = fused_ves(1:rxn_count);
                warning('SSA:ExceededCapacity',...
                        'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
                return;
            end

            T_cur = T_cur + tau_cur + tau;
            if collect_states
                T(rxn_count+1)   = T_cur; 
            end
            %%%%%%Carry out reaction
            %Note that the replenishment rate is not taken into account
            %here, as ves_state==pool_state yields that this is the only
            %possible reaction.
            
            if ves_ind > length(X)
                ves_act_ind = ves_ind - length(X);
                act_number_ves = act_number_vesicles(ves_act_ind);
                
                if abs(act_number_ves) == 1
                    act_number_vesicles(ves_act_ind) = 0;
                elseif act_number_ves == 0
                    
                
                [~, inact_rate] = determine_activation_rates(act_rate_const, inact_rate_const, Ca_R_temp, act_model_type); %alpha and beta
               
                    
                a_act_ves = [inact_rate, 0, delay_rate];
                react_number = find((cumsum(a_act_ves) >= r(3)*sum(a_act_ves)),1,'first'); % Determine reaction number;
                
                act_number_vesicles(ves_act_ind) = react_number - 2;
                end
                
                react_ind = 101; %Arbitrary number, set to avoid react_ind = 0 and thereby counting fusion twice
                
            elseif X(ves_ind) >= fuse_state
                a_site = propensities_newAlg_site_SS(X(ves_ind), par, Ca_R_temp(ves_ind), SS_PM, Ca_prim_type);
                react_number = find((cumsum(a_site) >= r(3)*sum(a_site)),1,'first'); % Determine reaction number;
                react_ind = react_number - 2;
                X(ves_ind) = X(ves_ind) + react_ind*10  + (abs(react_ind)-1)*100;
                react_ind = 101; %Arbitrary number, set to avoid react_ind = 0 and thereby counting fusion twice
            else 
                act_status = (act_number_vesicles(mod(ves_ind-1, num_ves)+1) == 1);
                [a_ves, n_ves, m_ves]  = propensities_newAlg_ves_SS(X(ves_ind), par, Ca_R_temp(ves_ind), act_status, Ca_prim_type); 

                react_number = find((cumsum(a_ves) >= r(3)*sum(a_ves)),1,'first'); % Determine reaction number;
            
                react_ind = react_number-3;
                
                X(ves_ind) = X(ves_ind)+(react_ind~=3)*(sign(react_ind).*10.^(abs(react_ind)-1)+(react_ind == 0).*((fuse_state-X(ves_ind)) + SS_PM*floor((X(ves_ind)/10))*10)) + (react_ind == 3)*(100 - n_ves + (SS_PM - 1)*m_ves);
            
                if ~collect_states && react_ind == 0
                    T(rxn_count+1) = T_cur;
                end
            
            
            end
            
            fused_ves(rxn_count+1) = (react_ind == 0);
            
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%
            if collect_states
                ves_states(:,rxn_count+1) = X;
                act_states(:,rxn_count+1) = act_number_vesicles;
            end
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%            
            
            
            if (collect_states==1)
                rxn_count = rxn_count + 1;
            elseif ~collect_states && react_ind == 0
                rxn_count = rxn_count + 1;
            end
            
            tau_cur = 0;

           
        end
    else %If all R-vesicles have fused
        T(rxn_count+1)   = tspan(2); 
        fused_ves(rxn_count+1) = 0;
        
        
        
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%
        if collect_states
            ves_states(:,rxn_count+1) = X;
            act_states(:,rxn_count+1) = act_number_vesicles;
        end
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%        
        

        rxn_count = rxn_count + 1;
    end
end

% Record output
t = T(1:rxn_count);
fused_ves = fused_ves(1:length(t));


%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%    
    if collect_states
        ves_states = ves_states(:,1:rxn_count);
        act_states = act_states(:,1:rxn_count);
    end
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%




if t(end) > tspan(2)
    t(end) = tspan(2);
    fused_ves(end) = 0;
    
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%  
    if collect_states
        ves_states(:,end) = X;
        act_states(:,end) = act_number_vesicles;
    end
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%    
    
end
if t(end) < tspan(2) %Because of tau_cur, algorithm can terminate without t(end)>=tspan(2), if t(end)+tau_cur >= tspan(2).
    t(end+1) = tspan(2);
    fused_ves(end+1) = 0;
    
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%      
    if collect_states
        ves_states(:,end+1) = X;
        act_states(:,end+1) = act_number_vesicles;
    end
%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%

end


if ~collect_states        
    ves_states = NaN; 
    act_states = NaN;
end

end