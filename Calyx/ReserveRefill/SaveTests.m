%[p_release_docked, p_release_tether, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, C_3, C_7]
lb = [0, 0, 0.0001, 0.0001, .0001, .0001, 0, 0, 0, 0];
ub = [0.6, 0.8, 0.01, 0.01, 0.01, 0.01, 100, 0.01, 0.5, 0.5];

%resolution = [5, 5, 4, 4, 4, 4, 2, 2, 2, 2];
resolution = [6, 2, 3, 2, 2, 3, 4, 2, 3, 6];
d = (ub-lb)./resolution;

for a = 1:resolution(1)
    x(1) = lb(1)+d(1)*a;
   
    for b = 1:resolution(2)
        x(2) = lb(2)+d(2)*b;
    
        for c = 1:resolution(3)
            x(3) = lb(3)+d(3)*c;
            
            for e = 1:resolution(4)
                x(4) = lb(4)+d(4)*e;
            
                for f = 1:resolution(5)
                    x(5) = lb(5)+d(5)*f;
                    
                    for g = 1:resolution(6)
                        x(6) = lb(6)+d(6)*g;

                        for h = 1:resolution(7)
                            x(7) = lb(7)+d(7)*h;

                            for i = 1:resolution(8)
                                x(8) = lb(8)+d(8)*i;

                                for j = 1:resolution(9)
                                    x(9) = lb(9)+d(9)*j;

                                    for k = 1:resolution(10)
                                        x(10) = lb(10)+d(10)*k;
                                        
                                        cost = DockingInhib2StateFunc(x);
                                        if cost < 200
                                            dataMat = [dataMat; [x cost]];
                                        else
                                            dataMat = [dataMat; [x nan]];
                                        end


                                    end 
                                end 
                            end 
                        end 
                    end
                end 
            end 
        end     
    end
end

save('SaveTests.mat','dataMat');


scatter3(dataMat(:,2),dataMat(:,3),dataMat(:,9),dataMat(:,11))
xlim = [0,1];
ylim = [0 0.1];
zlim = [0 1];
view(-31,14)
xlabel('p_{release_{tethered}}')
ylabel('k_{docking}')
zlabel('C_3')
cb = colorbar;
cb.Label.String = 'Cost';