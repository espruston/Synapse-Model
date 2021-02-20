clear
k_on_3 = 1e5; %M^-1ms^-1
k_off_3 = 0.2; %ms^-1  Hui
k_on_7 = 0.157; %Knight
k_off_7 = 0.011; %Knight
Ca_half_3 = 7e-6;
Ca_half_7 = 1.5e-6;
delta_t = 0.1;

SS_3 = 0.0244; %amount of bound syt3 in 50nm Ca
%make 1hz syt traces
mat = matfile('CF1HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat));
Syt3(1) = SS_3;

for i = 1:length(mat(1,:))-1

    Syt3(i+1) = Syt3(i) + (mat(i)*k_on_3*(1-Syt3(i))-k_off_3*Syt3(i))*delta_t;

end

figure('Name','Syt3 Simulation','NumberTitle','off')
subplot(4,1,1)
plot(Syt3)
title('1 Hz Syt3')
ylabel('Bound Syt3')

save('CF1HzSyt3.mat','Syt3');

%make 10hz syt traces
mat = matfile('CF10HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat));
Syt3(1) = SS_3; %SS val in 50nM

for i = 1:length(mat(1,:))-1

    Syt3(i+1) = Syt3(i) + (mat(i)*k_on_3*(1-Syt3(i))-k_off_3*Syt3(i))*delta_t;

end
subplot(4,1,2)
plot(Syt3)
title('10 Hz Syt3')
ylabel('Bound Syt3')

save('CF10HzSyt3.mat','Syt3');

%make 20hz syt traces
mat = matfile('CF20HzCa.mat').Ca_sim;
Syt3 = zeros(1,length(mat));
Syt3(1) = SS_3; %SS val in 50nM

for i = 1:length(mat(1,:))-1

    Syt3(i+1) = Syt3(i) + (mat(i)*k_on_3*(1-Syt3(i))-k_off_3*Syt3(i))*delta_t;

end
subplot(4,1,3)
plot(Syt3)
title('20 Hz Syt3')
ylabel('Bound Syt3')

save('CF20HzSyt3.mat','Syt3');

%make 50hz w/rec syt traces
mat = matfile('CF50HzCa.mat').Ca_sim;
Syt3 = zeros(size(mat));
Syt3(:,1) = SS_3; %SS val in 50nM

for j = 1:size(mat,1)
    for i = 1:size(mat,2)-1

        Syt3(j,i+1) = Syt3(j,i) + (mat(j,i)*k_on_3*(1-Syt3(j,i))-k_off_3*Syt3(j,i))*delta_t;

    end
end
subplot(4,1,4)
plot(Syt3(1,:))
title('50 Hz Syt3')
ylabel('Bound Syt3')

save('CF50HzSyt3.mat','Syt3');

% Syt7 = zeros(1,length(mat(1,:)));
% Syt7Ca(:) = mat(:).^2.8./(mat(:).^2.8+Ca_half_7.^2.8);
% 
% for i = 1:length(mat(1,:))-1
% 
%     %Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
%     %Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;
%     Syt7(i+1) = Syt7(i) + (1-Syt7(i))*Syt7Ca(i)*k_on_7*delta_t - Syt7(i)*k_off_7*delta_t;
%     
% end
% 
% 
% save('CF1HzSyt3.mat','Syt3');
% save('CF1HzSyt7.mat','Syt7');
% 
% %make 10hz syt traces
% mat = matfile('CF10HzCa.mat').Ca_sim;
% Syt3 = zeros(1,length(mat(1,:)));
% Syt7 = zeros(1,length(mat(1,:)));
% Syt3(:,1) = 0.0244; %SS val in 50nM
% %Syt7(:,1) = 0.0476;
% Syt7Ca = mat(:).^2.8./(mat(:).^2.8+Ca_half_7.^2.8);
% 
% for i = 1:length(mat(1,:))-1
% 
% %     Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
% %     Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;
%     Syt7(i+1) = Syt7(i) + (1-Syt7(i))*Syt7Ca(i)*k_on_7*delta_t - Syt7(i)*k_off_7*delta_t;
% end
% 
% save('CF10HzSyt3.mat','Syt3');
% save('CF10HzSyt7.mat','Syt7');
% 
% %make 20hz syt traces
% mat = matfile('CF20HzCa.mat').Ca_sim;
% Syt3 = zeros(1,length(mat(1,:)));
% Syt7 = zeros(1,length(mat(1,:)));
% Syt3(:,1) = 0.0244; %SS val in 50nM
% Syt7(:,1) = 0.0476;
% 
% for i = 1:length(mat(1,:))-1
% 
%     Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
%     Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;
% 
% end
% 
% save('CF20HzSyt3.mat','Syt3');
% save('CF20HzSyt7.mat','Syt7');
% 
% %make 50hz w/rec syt traces
% mat = matfile('CF50HzCa.mat').Ca_sim;
% Syt3 = zeros(1,length(mat(1,:)));
% Syt7 = zeros(1,length(mat(1,:)));
% Syt3(:,1) = 0.0244; %SS val in 50nM
% Syt7(:,1) = 0.0476;
% 
% for i = 1:length(mat(1,:))-1
% 
%     Syt3(:,i+1) = Syt3(:,i) + (mat(:,i).*k_on_3.*(1-Syt3(:,i))-k_off_3.*Syt3(:,i)).*delta_t;
%     Syt7(:,i+1) = Syt7(:,i) + (mat(:,i).*k_on_7.*(1-Syt7(:,i))-k_off_7.*Syt7(:,i)).*delta_t;
% 
% end
% 
% save('CF50HzSyt3.mat','Syt3');
% save('CF50HzSyt7.mat','Syt7');