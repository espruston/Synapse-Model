function result_mat = Two_Pool_Model_func(frequency, p_release_1, p_release_2, p_release_1_f, p_release_2_f, n_1_sites, n_2_sites, n_pulses, cpv, T_refill_1, T_refill_2, T_facilitation)

    delta_t = 1/frequency;
    n_1 = n_1_sites;
    n_2 = n_2_sites;
    
    n_1_empty = 0;
    n_2_empty = 0;

    n_1_f = 0;
    n_2_f = 0;
    
    n_1_new = 0;
    n_1_f_new = 0;
    n_2_new = 0;
    n_2_f_new = 0;
    
    n_released_1 = 0;
    n_released_2 = 0;
    n_released_1_f = 0;
    n_released_2_f = 0;
    
    n_released = 0;
    current = 0;
    
    result_mat = zeros(n_pulses + 1, 17);
    
    for i = 1:n_pulses
        
        result_mat(i,:) = [i-1, n_1_empty, n_2_empty, n_1, n_2, n_1_f, n_2_f, n_released_1, n_released_2, n_released_1_f, n_released_2_f, n_released, current, n_1_new, n_2_new, n_1_f_new, n_2_f_new]; 
        
        n_released_1 = n_1*p_release_1;
        n_released_2 = n_2*p_release_2;
        n_released_1_f = n_1_f*p_release_1_f;
        n_released_2_f = n_2_f*p_release_2_f;
        
        n_released = n_released_1 + n_released_2 + n_released_1_f + n_released_2_f;
        
        current = n_released*cpv;
        
        n_1_new = (n_1_empty + n_released_1 + n_released_1_f)*(1-exp(-1*delta_t/T_refill_1));
        n_2_new = (n_2_empty + n_released_2 + n_released_2_f)*(1-exp(-1*delta_t/T_refill_2));
        n_1_f_new = 0;%(n_1 - n_released_1)*(1-exp(-1*delta_t/T_facilitation));
        n_2_f_new = (n_2 - n_released_2)*(1-exp(-1*delta_t/T_facilitation));
        
        n_1 = n_1 - n_released_1 + n_1_new - n_1_f_new;
        n_1_f = n_1_f - n_released_1_f + n_1_f_new;
        
        n_2 = n_2 - n_released_2 + n_2_new - n_2_f_new;
        n_2_f = n_2_f - n_released_2_f + n_2_f_new;
        
        n_1_empty = n_1_sites - n_1 - n_1_f;
        n_2_empty = n_2_sites - n_2 - n_2_f;
        
    end

    result_mat(n_pulses + 1, :) = [n_pulses, n_1_empty, n_2_empty, n_1, n_2, n_1_f, n_2_f, n_released_1, n_released_2, n_released_1_f, n_released_2_f, n_released, current, n_1_new, n_2_new, n_1_f_new, n_2_f_new];

end