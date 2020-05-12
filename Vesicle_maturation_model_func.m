%Vesicle Maturation Model Function

function result_mat = Vesicle_maturation_model_func(v)
    % Vesicle maturation model, s1 = "lightly bound" vesicle, s2 = "strongly
    % bound" vesicle, s3 = ca facilitated vesicle
    frequency = v(1);
    T_refill = v(2);
    T_1_2 = v(3);
    T_facilitation = v(4);
    p_release_1 = v(5);
    p_release_2 = v(6);
    p_release_3 = v(7);
    n_vesicles = v(8);
    n_pulses = v(9);
    cpv = v(10);
    
    %test values
    %frequency = 100; %define pulse frequency (Hz)
    %T_1_2 = 7.5; %characteristic time (s) of maturation from stage 1->2 
    %T_refill = .25; %characteristic time (s) of t1 joining RRP
    %T_facilitation_2 = .025; %characteristic time (s) of t2 being facilitated by calcium influx due to pulse
    %p_release_1 = .058; %define probability of a "stage one" vesicle releasing
    %p_release_2 = .26; %define probability of a "stage two" vesicle releasing
    %p_release_3 = .4; %define probabiility of a "stage three" vesicle releasing
    %n_vesicles = 100; %define number of vesicles modeled
    %n_pulses = 50; %define number of pulses

    delta_t = 1/frequency; %time between pulses (s)

    n_s0 = 0; %no available docking sites-all sites have a vesicle to start
    n_s1 = n_vesicles; %all vesicles t1 to start before infinite time
    n_s2 = 0; %all vesicles t1 to start
    n_s3 = 0; %all vesicles t1 to start

    n_s0_new = 0; %new empty docking sites
    n_s1_new = 0; %number of t1 vesicles added after 'pulse 0'
    n_s2_new = n_vesicles; %number of t2 vesicles added after 'pulse 0' (infinite time)
    n_s3_new = 0; %number of t3 vesicles added after 'pulse 0'

    n_released_1 = 0; %number of t1 vesicles released
    n_released_2 = 0; %number of t2 vesicles released
    n_released_3 = 0; %number of t3 vesicles released

    n_released = 0; %total number of released vesicles is the sum of the number of each stage released

    current = 0; %abstract current per vesicle (stage independent)

    result_mat = zeros(n_pulses + 1, 14); %matrix with values per pulse

    for i = 1:n_pulses

        n_s1 = n_s1 - n_released_1 - n_s2_new + n_s1_new; %new total number of t1 vesicles at pulse i+1 = number at i - number released - number matured to stage 2 + number of newly bound
        n_s2 = n_s2 - n_released_2 - n_s3_new + n_s2_new; %new total number of t2 vesicles at pulse i+1
        n_s3 = n_s3 - n_released_3 +  + n_s3_new; %new total number of t3 vesicles at pulse i+1
        n_s0_new = n_released - n_s1_new;
        n_s0 = n_s0 + n_s0_new; %n_vesicles - n_s1 - n_s2 - n_s3; %new number of available docking sites

        result_mat(i,:) = [i - 1, n_s0, n_s1, n_s2, n_s3, n_released_1, n_released_2, n_released_3, n_released, current, n_s0_new, n_s1_new, n_s2_new, n_s3_new];

        n_released_1 = n_s1*p_release_1; %number of t1 vesicles released
        n_released_2 = n_s2*p_release_2; %number of t2 vesicles released
        n_released_3 = n_s3*p_release_3; %number of t3 vesicles released

        n_released = n_released_1 + n_released_2 + n_released_3; %total number of released vesicles is the sum of the number of each stage released

        current = n_released*cpv; %abstract current of unit -1 per vesicle (stage independent)

        n_s1_new = (n_s0 + n_released)*(1 - exp(-1*delta_t/T_refill)); %the number of t1 vesicles which join the pool after release from pulse n and before pulse n+1 is dependent upon the number of available docking sites and the logarithmic growth of the pool in delta_t seconds with characteristic time T_refill

        n_s2_new = (n_s1 - n_released_1)*(1 - exp(-1*delta_t/T_1_2)); %the number of t1 vesicles that become t2 vesicles after pulse n and before pulse n+1 is modeled as logarithmic growth with characteristic time T_1_2 in time delta_t

        n_s3_new = (n_s2 - n_released_2)*(1 - exp(-1*delta_t/T_facilitation)); %0 represents no facilitation of stage 1 vesicles; the number of t2 vesicles that become facilitated after pulse n and before pusle n+1 is modeled as logarithmic growth in delta_t time with characteristic time T_facilitation


    end
    result_mat(n_pulses+1,:) = [n_pulses, n_s0, n_s1, n_s2, n_s3, n_released_1, n_released_2, n_released_3, n_released, current, n_s0_new, n_s1_new, n_s2_new, n_s3_new];

end