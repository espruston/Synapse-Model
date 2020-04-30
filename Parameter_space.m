% parameter space probe script

%vesicle maturation

%SET THESE PARAMETERS FOR DESIRED PARAMETER SPACE
n_vesicles = 100;
n_pulses = 100;
cpv = -1;
frequency = 100;

t_1 = (0.01:0.01:12); %range of stage 0 -> 1 time constant (s) AKA T_REFILL
t_2 = (0.01:0.01:12); %range of stage 1 -> 2 time constant (s) AKA T_MATURATION
t_3 = (0.01:0.01:12); %range of stage 2 -> 3 time constant (s) AKA T_FACILITATION

p_1 = (0.001:0.001:1); %range of probabilities for release of stage 1 AKA immature release
p_2 = (0.001:0.001:1); %range of probabilities for release of stage 2 AKA mature release
p_3 = (0.001:0.001:1); %range of probabilities for release of stage 3 AKA facilitated release

len_t_1 = size(t_1,2);
len_t_2 = size(t_2,2);
len_t_3 = size(t_3,2);

len_p_1 = size(p_1,2);
len_p_2 = size(p_2,2);
len_p_3 = size(p_3,2);

n_results = len_t_1 * len_t_2 * len_t_3 * len_p_1 * len_p_2 * len_p_3;

promptMessage = sprintf('WARNING: THIS PROGRAM CONTAINS A LOOP THAT RUNS OVER %d ITERATIONS, ARE YOU SURE YOU WANT TO CONTINUE?', n_results);
button = questdlg(promptMessage, 'Continue', 'Continue', 'Abort', 'Continue');
if strcmp(button, 'Abort')
  return; % or break or whatever...
end

for i1 = 1:len_t_1
    T_refill = t_1(i1);
    
    for i2 = 1:len_t_2
        T_1_2 = t_2(i2);
        
        for i3 = 1:len_t_3
            T_facilitation = t_3(i3);
            
            for j1 = 1:len_p_1
                p_release_1 = p_1(j1);
                
                for j2 = 1:len_p_2
                    p_release_2 = p_2(j2);
                    
                    for j3 = 1:len_p_3
                        p_release_3 = p_3(j3);
                        
                        
                        Vesicle_maturation_model_func(frequency, T_refill, T_1_2, T_facilitation, p_release_1, p_release_2, p_release_3, n_vesicles, n_pulses, cpv)
    
                        
    
    
    
    
                    end
                end
            end
        end
    end
end
