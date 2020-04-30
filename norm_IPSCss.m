two_pool_current = zeros(10,101);
vesicle_maturation_current = zeros(10,101);
Hz = zeros(1,10);
for i = 1:10
    
    f = i*10;
    
    two_pool_mat = Two_Pool_Model_func(f, .14, 0.023, 0, 0.035, 30, 70, 100, -1, 7.5, 0.28, 0.125);
    two_pool_current(i,:) = two_pool_mat(:,13)/min(two_pool_mat(:,13)); %normalized current
    
    vesicle_maturation_mat = Vesicle_maturation_model_func(f, .2, 7.5, 0.025, .028, .09, .09, 100, 100, -1);
    vesicle_maturation_current(i,:) = vesicle_maturation_mat(:,10)/min(vesicle_maturation_mat(:,10));
    
    Hz(i) = f;
end

F1 = line(Hz,two_pool_current(:,101));
F1.Color = 'red';
F1.LineStyle = '-';
F1.Marker = 'o';
title("Frequency dependence of IPSC_{ss}")
xlabel("Frequency (Hz)")
ylabel("IPSC_{ss} (norm)")
ylim([0,.6])

hold on
grid on

F2 = line(Hz,vesicle_maturation_current(:,101));
F2.Color = 'blue';
F2.LineStyle = '-';
F2.Marker = 'x';

legend('Two Pool', 'Vesicle Maturation')
