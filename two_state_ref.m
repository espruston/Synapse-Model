%two state vectorized ODE reference file

k_AB = 1;
k_BA = 2;

rate_mat = [0, k_BA; k_AB, 0];

rate_mat = rate_mat + -1*eye(2).*sum(rate_mat); %create rate matrix where the on diagonal components are equal to -1 times the sum of all outward rate constants for A

%this method extends to an infinite # of states where rate_mat_i,j = k_ij

ICs = [10; 20]; %initial conditions, [A; B]
time_interval = [0,100];

written_wrapper = @(t,y) ODEs_written(t,y,k_AB,k_BA); %wrapper function for passing k values to ode solver
[t_written, y_written] = ode45(written_wrapper, time_interval, ICs);
vec_wrapper = @(t,y) ODEs_vec(t,y,rate_mat); %wrapper function for passing rate matrix to ode solver
[t_vec, y_vec] = ode45(vec_wrapper, time_interval, ICs);

subplot(2,1,1);
plot(t_written, y_written(:,1),t_written,y_written(:,2))
legend('A','B')
title('Written solution');
xlabel('Time');
ylabel('y');

subplot(2,1,2);
plot(t_vec, y_vec(:,1),t_vec,y_vec(:,2))
legend('A','B')
title('Vectorized solution');
xlabel('Time');
ylabel('y');

function dydt = ODEs_written(t,y,k_AB,k_BA)
    
    %y = [A, B]^T
    dAdt = y(2)*k_BA - y(1)*k_AB;
    dBdt = y(1)*k_AB - y(2)*k_BA;
    
    dydt = [dAdt; dBdt];

end

function dydt = ODEs_vec(t,y,rate_mat)
    
    %y = [A,B]^T
   
    dydt = [rate_mat*y];
    
end


