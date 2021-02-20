clear
k_on_3 = 1e5; %M^-1ms^-1
k_off_3 = 0.2; %ms^-1  Hui
k_on_7 = 1e4; %M^-1ms^-1
k_off_7 = 0.010; %Knight
Ca_half_3 = 7e-6;
Ca_half_7 = 1.5e-6;
delta_t = 0.1;

SS_3 = 0.0244; %amount of bound syt3 in 50nm Ca
SS_7 = 0.0476;

%make 1hz syt traces
mat = matfile('Calyx1HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat));
Syt7 = zeros(1,length(mat));
Syt3(1) = SS_3;

for i = 1:length(mat(1,:))-1

    Syt3(i+1) = Syt3(i) + (mat(i)*k_on_3*(1-Syt3(i))-k_off_3*Syt3(i))*delta_t;
    Syt7(i+1) = Syt7(i) + (mat(i)*k_on_7*(1-Syt7(i))-k_off_7*Syt7(i))*delta_t;

end

figure('Name','Syt Simulation','NumberTitle','off')
subplot(6,1,1)
plot(Syt3)
hold on
plot(Syt7)
title('1 Hz Syt')
ylabel('Bound Syt')

save('Calyx1HzSyt3.mat','Syt3');
save('Calyx1HzSyt7.mat','Syt7');

%make 10hz w/rec syt traces
mat = matfile('Calyx10HzCa.mat').Ca_sim;
Syt3 = zeros(size(mat));
Syt7 = zeros(size(mat));
Syt3(:,1) = SS_3; %SS val in 50nM
Syt7(:,1) = SS_7; %SS val in 50nM

for j = 1:size(mat,1)
    for i = 1:size(mat,2)-1

        Syt3(j,i+1) = Syt3(j,i) + (mat(j,i)*k_on_3*(1-Syt3(j,i))-k_off_3*Syt3(j,i))*delta_t;
        Syt7(j,i+1) = Syt7(j,i) + (mat(j,i)*k_on_7*(1-Syt7(j,i))-k_off_7*Syt7(j,i))*delta_t;

    end
end
subplot(6,1,2)
plot(Syt3(1,:))
hold on
plot(Syt7(1,:))
title('10 Hz Syt3')
ylabel('Bound Syt3')

save('Calyx10HzSyt3.mat','Syt3');
save('Calyx10HzSyt7.mat','Syt7');

%make 20hz w/rec syt traces
mat = matfile('Calyx20HzCa.mat').Ca_sim;
Syt3 = zeros(size(mat));
Syt7 = zeros(size(mat));
Syt3(:,1) = SS_3; %SS val in 50nM
Syt7(:,1) = SS_7; %SS val in 50nM

for j = 1:size(mat,1)
    for i = 1:size(mat,2)-1

        Syt3(j,i+1) = Syt3(j,i) + (mat(j,i)*k_on_3*(1-Syt3(j,i))-k_off_3*Syt3(j,i))*delta_t;
        Syt7(j,i+1) = Syt7(j,i) + (mat(j,i)*k_on_7*(1-Syt7(j,i))-k_off_7*Syt7(j,i))*delta_t;

    end
end
subplot(6,1,3)
plot(Syt3(1,:))
hold on
plot(Syt7(1,:))
title('20 Hz Syt')
ylabel('Bound Syt')

save('Calyx20HzSyt3.mat','Syt3');
save('Calyx20HzSyt7.mat','Syt7');

%make 50hz w/rec syt traces
mat = matfile('Calyx50HzCa.mat').Ca_sim;
Syt3 = zeros(size(mat));
Syt7 = zeros(size(mat));
Syt3(:,1) = SS_3; %SS val in 50nM
Syt7(:,1) = SS_7; %SS val in 50nM

for j = 1:size(mat,1)
    for i = 1:size(mat,2)-1

        Syt3(j,i+1) = Syt3(j,i) + (mat(j,i)*k_on_3*(1-Syt3(j,i))-k_off_3*Syt3(j,i))*delta_t;
        Syt7(j,i+1) = Syt7(j,i) + (mat(j,i)*k_on_7*(1-Syt7(j,i))-k_off_7*Syt7(j,i))*delta_t;

    end
end
subplot(6,1,4)
plot(Syt3(1,:))
hold on
plot(Syt7(1,:))
title('50 Hz Syt')
ylabel('Bound Syt')

save('Calyx50HzSyt3.mat','Syt3');
save('Calyx50HzSyt7.mat','Syt7');

%make 100hz w/rec syt traces
mat = matfile('Calyx100HzCa.mat').Ca_sim;
Syt3 = zeros(size(mat));
Syt7 = zeros(size(mat));
Syt3(:,1) = SS_3; %SS val in 50nM
Syt7(:,1) = SS_7; %SS val in 50nM

for j = 1:size(mat,1)
    for i = 1:size(mat,2)-1

        Syt3(j,i+1) = Syt3(j,i) + (mat(j,i)*k_on_3*(1-Syt3(j,i))-k_off_3*Syt3(j,i))*delta_t;
        Syt7(j,i+1) = Syt7(j,i) + (mat(j,i)*k_on_7*(1-Syt7(j,i))-k_off_7*Syt7(j,i))*delta_t;

    end
end
subplot(6,1,5)
plot(Syt3(1,:))
hold on
plot(Syt7(1,:))
title('100 Hz Syt')
ylabel('Bound Syt')

save('Calyx100HzSyt3.mat','Syt3');
save('Calyx100HzSyt7.mat','Syt7');

%make 200hz w/rec syt traces
mat = matfile('Calyx200HzCa.mat').Ca_sim;
Syt3 = zeros(size(mat));
Syt7 = zeros(size(mat));
Syt3(:,1) = SS_3; %SS val in 50nM
Syt7(:,1) = SS_7; %SS val in 50nM

for j = 1:size(mat,1)
    for i = 1:size(mat,2)-1

        Syt3(j,i+1) = Syt3(j,i) + (mat(j,i)*k_on_3*(1-Syt3(j,i))-k_off_3*Syt3(j,i))*delta_t;
        Syt7(j,i+1) = Syt7(j,i) + (mat(j,i)*k_on_7*(1-Syt7(j,i))-k_off_7*Syt7(j,i))*delta_t;

    end
end
subplot(6,1,6)
plot(Syt3(1,:))
hold on
plot(Syt7(1,:))
title('200 Hz Syt')
ylabel('Bound Syt')

save('Calyx200HzSyt3.mat','Syt3');
save('Calyx200HzSyt7.mat','Syt7');
