function [steadyvesicles_reshaped] = calculate_steady_state(par, Ca_bas, num_ves_in_steady)

%Determines the steady state of calcium binding to one or two calcium
%sensors. Steady state of activation or priming is determined in
%determining_states_with_fusion.m


k_4 = par(6);
k_min4 = par(7);
b_s = par(8);
k_3 = par(1);
k_min3 = par(2);
b_f = par(3);

n_max = par(13);
m_max = par(14);
SS_coop = m_max;


steady_state = ones(m_max+1, n_max+1);

if SS_coop > 0  %These are models with slow sensor
    for n = 0:n_max
        for m = 0: m_max

            if m == 0 && n > 0
                steady_state(1,n+1) = ((factorial(n_max)/factorial(n_max-n)) * Ca_bas^n * k_3^n)/(factorial(n)*(b_f^(n*(n-1)/2))*k_min3^n);      
            elseif n == 0 && m > 0
                steady_state(m+1,1) = ((factorial(m_max)/factorial(m_max-m)) * Ca_bas^m * k_4^m)/(factorial(m)*(b_s^(m*(m-1)/2))*k_min4^m);         
            elseif n > 0 && m > 0
                steady_state(m+1, n+1) = (((factorial(n_max)/factorial(n_max-n)) * Ca_bas^n * k_3^n)/(factorial(n)*(b_f^(n*(n-1)/2))*k_min3^n))  *  (((factorial(m_max)/factorial(m_max-m)) * Ca_bas^m * k_4^m)/(factorial(m)*(b_s^(m*(m-1)/2))*k_min4^m));  
            end
        end
    end
elseif SS_coop == 0 %These are models without slow sensor
    for n = 0:n_max
        for m = 0: m_max

            if m == 0 && n > 0
                steady_state(1,n+1) = ((factorial(n_max)/factorial(n_max-n)) * Ca_bas^n * k_3^n)/(factorial(n)*(b_f^(n*(n-1)/2))*k_min3^n);      
            elseif n == 0 && m > 0
                steady_state(m+1,1) = 0;         
            elseif n > 0 && m > 0
                steady_state(m+1, n+1) = 0;  
            end
        end
    end
else
    disp('WARNING. SS_on_off VALUE IS INVALID')
    return
end
    
    
    
 
steadyvesicles_reshaped_first = reshape(steady_state,[],1);
steady_factor = num_ves_in_steady/sum(steadyvesicles_reshaped_first);
steadyvesicles_reshaped = (steadyvesicles_reshaped_first*steady_factor);
    
    
    
    
    
  