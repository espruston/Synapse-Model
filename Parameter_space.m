% parameter space probe script

%vesicle maturation

%SET THESE PARAMETERS FOR DESIRED PARAMETER SPACE
n_vesicles = 100;
n_pulses = 100;
cpv = -1;

frequency = 100; %range of frequencies probed

t_1 = (0.01:0.01:12); %range of stage 0 -> 1 time constant (s) AKA T_REFILL
t_2 = (0.01:0.01:12); %range of stage 1 -> 2 time constant (s) AKA T_MATURATION
t_3 = (0.01:0.01:12); %range of stage 2 -> 3 time constant (s) AKA T_FACILITATION

p_1 = (0.001:0.001:1); %range of probabilities for release of stage 1 AKA immature release
p_2 = (0.001:0.001:1); %range of probabilities for release of stage 2 AKA mature release
p_3 = (0.001:0.001:1); %range of probabilities for release of stage 3 AKA facilitated release

len_f = size(frequencies,2);

len_t_1 = size(t_1,2);
len_t_2 = size(t_2,2);
len_t_3 = size(t_3,2);

len_p_1 = size(p_1,2);
len_p_2 = size(p_2,2);
len_p_3 = size(p_3,2);