function [SS_propensities] = calculate_SS_propensities(Ca_R_vesicles, input_population, par, Ca_time, t)

%%%%THIS SCRIPT GENERATES PROPENSITIES FOR THE CALCIUM BINDING TO THE
%%%%SECOND SENSOR REGARDLESS OF THE FIRST SENSOR. IS USED FOR CALCIUM
%%%%BINDING IN EMPTY SITES AND CAN BE USED FOR POPULATED STATES.


Calcium = interp1(Ca_time, Ca_R_vesicles, t);

SS_propensities = zeros(size(input_population));

input_size = size(input_population);
states_length = input_size(2);


%%%%PARAMETERS

% k_3  = par(1);
% k_3b = par(2);
% b_f  = par(3);
% l_0 = par(4);
% f = par(5);
k_4  = par(6);
k_4b = par(7);
b_s  = par(8);
% s = par(9);
n_max = par(13);
m_max = par(14);

% if f == 0 && s == 0
%     l_0 = 0;
% end
% 
% if allow_fus == 0
%     f = 1;
%     s = 1;
%     l_0 = 0;
% end



%%%%





if m_max > 0

    for m = 0:m_max
        for n = 0:(states_length-1)

            if m > 0 && m < m_max && n > 0 && n < n_max
                % for 0<m<mmax && 0<n<nmax % general case ('middle middle')
                SS_propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1);       %[S(m)F(n)]


            
            elseif m == 0 && n == 0
                % for m = 0 & n = 0 % top left
                 SS_propensities(m+1,n+1) = (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1);       %[S(m)F(n)]  

                

            elseif m == 0 && n > 0 && n < n_max
                % for m = 0, 0 < n < n_max % top center
                SS_propensities(m+1,n+1) =  (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1);       %[S(m)F(n)]  



            elseif m == 0 && n == n_max
                % for m = 0, n = n_max % top right
                            SS_propensities(m+1,n+1) =   (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1);       %[S(m)F(n)]

 
            elseif m > 0 && m < m_max && n == 0
                       SS_propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1);       %[S(m)F(n)]

 
            elseif m > 0 && m < m_max && n == n_max
                % for 0<m<m_max, n = n_max % middle right
                          SS_propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                + (m+1)         * k_4b * b_s^m      * input_population(m+1+1,n+1) ...  %[S(m+1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1) ...    %[S(m)F(n)]
                                - (m_max - m)   * k_4 * Calcium     * input_population(m+1,n+1);       %[S(m)F(n)]


            elseif m == m_max && n == 0
                % for m = m_max, n = 0 % bottom left
                            SS_propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1);    %[S(m)F(n)]

 
            elseif m == m_max && n > 0 && n < n_max
                % for m = m_max, 0<n<n_max % bottom middle
                            SS_propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1);    %[S(m)F(n)]


  
            elseif m == m_max && n == n_max
                %for m = m_max, n = n_max % bottom right
    %              dual_SS(m+1,n+1) =    num_ves - sum(sum(input_pop))-input_pop(m+1,n+1); % Mass conservation
                          SS_propensities(m+1,n+1) =    (m_max - m+1) * k_4 * Calcium     * input_population(m+1-1,n+1) ...  %[S(m-1)F(n)]
                                 - m             * k_4b * b_s^(m-1)  * input_population(m+1,n+1);    %[S(m)F(n)]

            else

                disp('ERROR - state not assigned');

            end
        end
    end
end