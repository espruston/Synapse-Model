clear
k_on_3 = 1e5; %M^-1ms^-1 Hui
k_off_3 = 0.2; %ms^-1  Hui
k_on_7 = 1e4; %Knight
k_off_7 = 0.01; %Knight
delta_t = 0.01;

%make 1hz syt traces
mat = matfile('CF1HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat(1,:)));
Syt7 = zeros(1,length(mat(1,:)));
Syt3(:,1) = 0.0244; %SS val in 50nM
Syt7(:,1) = 0.0476;

for i = 1:length(mat(1,:))-1

    Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
    Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;

end

save('CF1HzSyt3.mat','Syt3');
save('CF1HzSyt7.mat','Syt7');

%make 10hz syt traces
mat = matfile('CF10HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat(1,:)));
Syt7 = zeros(1,length(mat(1,:)));
Syt3(:,1) = 0.0244; %SS val in 50nM
Syt7(:,1) = 0.0476;

for i = 1:length(mat(1,:))-1

    Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
    Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;

end

save('CF10HzSyt3.mat','Syt3');
save('CF10HzSyt7.mat','Syt7');

%make 20hz syt traces
mat = matfile('CF20HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat(1,:)));
Syt7 = zeros(1,length(mat(1,:)));
Syt3(:,1) = 0.0244; %SS val in 50nM
Syt7(:,1) = 0.0476;

for i = 1:length(mat(1,:))-1

    Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
    Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;

end

save('CF20HzSyt3.mat','Syt3');
save('CF20HzSyt7.mat','Syt7');

%make 50hz w/rec syt traces
mat = matfile('CF50HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat(1,:)));
Syt7 = zeros(1,length(mat(1,:)));
Syt3(:,1) = 0.0244; %SS val in 50nM
Syt7(:,1) = 0.0476;

for i = 1:length(mat(1,:))-1

    Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
    Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;

end

save('CF50HzSyt3.mat','Syt3');
save('CF50HzSyt7.mat','Syt7');