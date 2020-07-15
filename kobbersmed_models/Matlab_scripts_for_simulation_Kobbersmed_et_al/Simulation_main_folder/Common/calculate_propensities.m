function [propensities, fusion_rate] = calculate_propensities(Ca_R_vesicles, input_population, par, Ca_time, t, allow_fus, SS_PM)

%%%%THIS SCRIPT GENERATES PROPENSITIES FOR THE CALCIUM BINDING AND FUSION
%%%%OF VESICLES

%Parameters
k_3  = par(1);
k_3b = par(2);
b_f  = par(3);
l_0 = par(4);
f = par(5);
k_4  = par(6);
k_4b = par(7);
b_s  = par(8);
s = par(9);
n_max = par(13);
m_max = par(14);

if f == 0 || s == 0
    l_0 = 0;
end

if allow_fus == 0
    f = 1;
    s = 1;
    l_0 = 0;
end

Calcium = interp1(Ca_time, Ca_R_vesicles, t);

propensities = zeros(size(input_population));

fusion_rate = 0;
if SS_PM == 1
    fusion_rate = zeros(m_max + 1, 1);
end


%%%%





if m_max > 0

    for m = 0:m_max
        for n = 0:n_max

            if m > 0 && m < m_max && n > 0 && n < n_max
                % for 0<m<mmax && 0<n<nmax % general case ('middle middle')
                propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (n_max - n+1) * k_3 * Calcium     * input_population(m+1,n+1-1) ...  %[S(m)F(n-1)]
                                + (n+1)         * k_3b * b_f^n      * input_population(m+1,n+1+1) ...  %[S(m)F(n+1)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - n             * k_3b * b_f^(n-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);       %[S(m)F(n)]
                                
                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end
            
            elseif m == 0 && n == 0
                % for m = 0 & n = 0 % top left
                 propensities(m+1,n+1) =     (n+1)         * k_3b * b_f^n      * input_population(m+1,n+1+1) ...  %[S(m)F(n+1)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);       %[S(m)F(n)]  

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end
                

            elseif m == 0 && n > 0 && n < n_max
                % for m = 0, 0 < n < n_max % top center
                propensities(m+1,n+1) =  (n_max - n+1) * k_3 * Calcium     * input_population(m+1,n+1-1) ...  %[S(m)F(n-1)]
                                + (n+1)         * k_3b * b_f^n      * input_population(m+1,n+1+1) ...  %[S(m)F(n+1)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - n             * k_3b * b_f^(n-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);       %[S(m)F(n)]  

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end

            elseif m == 0 && n == n_max
                % for m = 0, n = n_max % top right
                            propensities(m+1,n+1) =   (n_max - n+1) * k_3 * Calcium     * input_population(m+1,n+1-1) ...  %[S(m)F(n-1)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - n             * k_3b * b_f^(n-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);       %[S(m)F(n)]

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end

 
            elseif m > 0 && m < m_max && n == 0
                       propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (n+1)         * k_3b * b_f^n      * input_population(m+1,n+1+1) ...  %[S(m)F(n+1)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);       %[S(m)F(n)]

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end
 
            elseif m > 0 && m < m_max && n == n_max
                % for 0<m<m_max, n = n_max % middle right
                          propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (n_max - n+1) * k_3 * Calcium     * input_population(m+1,n+1-1) ...  %[S(m)F(n-1)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - n             * k_3b * b_f^(n-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);       %[S(m)F(n)]

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end


            elseif m == m_max && n == 0
                % for m = m_max, n = 0 % bottom left
                            propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (n+1)         * k_3b * b_f^n      * input_population(m+1,n+1+1) ...  %[S(m)F(n+1)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);    %[S(m)F(n)]

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end

 
            elseif m == m_max && n > 0 && n < n_max
                % for m = m_max, 0<n<n_max % bottom middle
                            propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (n_max - n+1) * k_3 * Calcium     * input_population(m+1,n+1-1) ...  %[S(m)F(n-1)]
                                + (n+1)         * k_3b * b_f^n      * input_population(m+1,n+1+1) ...  %[S(m)F(n+1)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - n             * k_3b * b_f^(n-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);    %[S(m)F(n)]

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end
  
            elseif m == m_max && n == n_max
                %for m = m_max, n = n_max % bottom right
    %              dual_SS(m+1,n+1) =    num_ves - sum(sum(input_pop))-input_pop(m+1,n+1); % Mass conservation
                          propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                 + (n_max - n+1) * k_3 * Calcium     * input_population(m+1,n+1-1) ...  %[S(m)F(n-1)]
                                 - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                 - n             * k_3b * b_f^(n-1)  * input_population(m+1,n+1) ...       %[S(m)F(n)]
                                - l_0 * s^m * f^n * input_population(m+1,n+1);    %[S(m)F(n)]

                            if SS_PM == 0
                                fusion_rate = fusion_rate + l_0 * s^m * f^n * input_population(m+1,n+1);
                            elseif SS_PM == 1
                                fusion_rate(m+1) = fusion_rate(m+1) + l_0 * s^m * f^n * input_population(m+1,n+1);
                            end
            else

                disp('ERROR - state not assigned');

            end
        end
    end
    
elseif m_max == 0 %IF SLOW SENSOR IS OFF
    for n = 0:n_max

            if n > 0 && n < n_max
                % for 0<m<mmax && 0<n<nmax % general case ('middle middle')
                propensities(n+1) = (n_max - n+1) * k_3 * Calcium     * input_population(n+1-1) ...  %[S(m)F(n-1)]
                                + (n+1)         * k_3b * b_f^n      * input_population(n+1+1) ...  %[S(m)F(n+1)]
                                - n             * k_3b * b_f^(n-1)  * input_population(n+1) ...    %[S(m)F(n)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(n+1) ...    %[S(m)F(n)]
                                - l_0 * f^n * input_population(n+1);       %[S(m)F(n)]

                                fusion_rate = fusion_rate + l_0 * f^n * input_population(n+1);
            elseif n == 0
                % for m = 0 & n = 0 % top left
                 propensities(n+1) =     (n+1)         * k_3b * b_f^n      * input_population(n+1+1) ...  %[S(m)F(n+1)]
                                - (n_max - n)   * k_3 * Calcium     * input_population(n+1) ...     %[S(m)F(n)]
                                - l_0 * f^n * input_population(n+1);       %[S(m)F(n)]  

                                fusion_rate = fusion_rate + l_0 * f^n * input_population(n+1);

                
            elseif n == n_max
                % for m = 0, n = n_max % top right
                            propensities(n+1) =   (n_max - n+1) * k_3 * Calcium     * input_population(n+1-1) ...  %[S(m)F(n-1)]
                                - n             * k_3b * b_f^(n-1)  * input_population(n+1) ...    %[S(m)F(n)]
                                - l_0 * f^n * input_population(n+1);       %[S(m)F(n)]

                                fusion_rate = fusion_rate + l_0 * f^n * input_population(n+1);

            end
    end
    
    
end