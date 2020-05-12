% parameter space probe script

%vesicle maturation
lambda = 9; %lambda value which is being matched

n_free_parameters = 6; %number of free parameters to be searhed
%SET THESE PARAMETERS FOR DESIRED PARAMETER SPACE
n_vesicles = 100;
n_pulses = 100;
cpv = 1;
frequency = 100;

t_1 = (0.01:0.1:3.01); %range of stage 0 -> 1 time constant (s) AKA T_REFILL
t_2 = (1:0.1:10.1); %range of stage 1 -> 2 time constant (s) AKA T_MATURATION
t_3 = (1:0.1:10.1); %range of stage 2 -> 3 time constant (s) AKA T_FACILITATION

p_1 = (0.001:0.005:.151); %range of probabilities for release of stage 1 AKA immature release
p_2 = (0.001:0.01:.251); %range of probabilities for release of stage 2 AKA mature release
p_3 = (0.001:0.01:.551); %range of probabilities for release of stage 3 AKA facilitated release

free_parameter_sizes = [size(t_1,2), size(t_2,2), size(t_3,2), size(p_1,2), size(p_2,2), size(p_3,2)]; %vector containing the size of each parameter range
parameter_mins = [min(t_1), min(t_2), min(t_3), min(p_1), min(p_2), min(p_3)];
parameter_maxs = [max(t_1), max(t_2), max(t_3), max(p_1), max(p_2), max(p_3)];

index_mat = zeros(n_free_parameters, max(free_parameter_sizes)); %create matrix of size free_parameters*maxsize

for i = 1:n_free_parameters
   
    n = free_parameter_sizes(i);
    index_mat(i,1:n) = randperm(n); %populate matrix with random indicies
    
end

post_GD_indicies = zeros(min(free_parameter_sizes), n_free_parameters+1); %create an emptpy matrix that will store the post search value for best indicies and variation from expected value
post_GD_values = zeros(min(free_parameter_sizes), n_free_parameters+1); %create an emptpy matrix that will store the post search value for best indicies and variation from expected value

for i = 1:min(free_parameter_sizes)
    
lowest_error = abs(lambda-n_pulses+1)+1; %one worse error than worst possible error
delta_lowest_error = 1;

indicies = index_mat(:,i);

    %run a gradient descent search for the values matching lambda
    while delta_lowest_error ~= 0 %as long as progress is being made

        lowest_error_start = lowest_error;

        tempvec = indicies-1; %subtract 1 from all indicies
        indicies_minus = tempvec + 1*(tempvec<1); %add one if the minimum index is passed

        tempvec = indicies+1;
        indicies_plus = tempvec - 1*(tempvec>free_parameter_sizes'); %subtract 1 from any index that exceeds bounds

        %create the matrix containing the parameter space
        param_space = [indicies_minus,indicies,indicies_plus]';
        param_space = combvec(param_space(:,1)',param_space(:,2)',param_space(:,3)',param_space(:,4)',param_space(:,5)',param_space(:,6)'); %n_free_parameters*3^n matrix of combinations of indicies

        for j=1:size(param_space,2) %iterate through all parameter space to find the best solution in the range

            input = [frequency, t_1(param_space(1,j)), t_2(param_space(2,j)), t_3(param_space(3,j)), p_1(param_space(4,j)), p_2(param_space(5,j)), p_3(param_space(6,j)), n_vesicles, n_pulses, cpv];
            temp_results = Vesicle_maturation_model_func(input);
            temp_currents = temp_results(:,10);
            current_drop = temp_currents(2) - temp_currents(n_pulses+1); %calculate the drop in current from pulse 1 to last pulse
            target_current = exp(-1)*current_drop + temp_currents(n_pulses+1); %the current value which is equal to 1/e*drop+steady state

            for k = 1:n_pulses

                if temp_currents(k+1) < target_current %compute the lambda value for the jth set of parameters
                    break
                end

            end

            error = abs(lambda - k+1);

            if error <= lowest_error
                lowest_error = error;
                j_lowest_error = j;
            end

        end

        delta_lowest_error = lowest_error_start - lowest_error;
        indicies = param_space(:,j_lowest_error); %placed at the end for easy recovery after loop
    end

post_GD_indicies(i,:) = [lowest_error, indicies'];
post_GD_values(i,:) = [lowest_error, t_1(indicies(1)), t_2(indicies(2)), t_3(indicies(3)), p_1(indicies(4)), p_2(indicies(5)), p_3(indicies(6))];
    
end