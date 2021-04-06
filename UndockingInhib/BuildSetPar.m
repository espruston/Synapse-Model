%[p_release_docked, p_release_tether, k_docking, k_undocking, k_tether, k_untether, reserve_size, k_refill, C_3, C_7]
lb = [0, 0, 0.0001, 0.0001, .0001, .0001, 0, 0, 0, 0]; %
ub = [0.6, 0.8, 0.01, 0.01, 0.01, 0.01, 100, 0.01, 0.5, 0.5];

resolution = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]; %each entry must be at least 2

size = 1;
for n = 1:length(resolution)
   
    size = size*resolution(n);
    
end

dM = zeros(size, 10); %matrix containing parameter set
dataMat = zeros(size, 1); %matrix containing parameter set and cost
traces = zeros(111,24,size); %matrix containing all 24 traces for each parameter set

m = 1; %dataMat and traces indexing

d = (ub-lb)./(resolution - 1);
x = lb;

%build params
for a = 1:resolution(1)
    x(1) = lb(1)+d(1)*(a-1);
   
    for b = 1:resolution(2)
        x(2) = lb(2)+d(2)*(b-1);
    
        for c = 1:resolution(3)
            x(3) = lb(3)+d(3)*(c-1);
            
            for e = 1:resolution(4)
                x(4) = lb(4)+d(4)*(e-1);
            
                for f = 1:resolution(5)
                    x(5) = lb(5)+d(5)*(f-1);
                    
                    for g = 1:resolution(6)
                        x(6) = lb(6)+d(6)*(g-1);

                        for h = 1:resolution(7)
                            x(7) = lb(7)+d(7)*(h-1);

                            for i = 1:resolution(8)
                                x(8) = lb(8)+d(8)*(i-1);

                                for j = 1:resolution(9)
                                    x(9) = lb(9)+d(9)*(j-1);

                                    for k = 1:resolution(10)
                                        x(10) = lb(10)+d(10)*(k-1);
                                        
                                        dM(m,:) = x;
                                        
                                        m = m+1;
                                        
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

%parallel test params
tic
parfor m = 1:size        
    
    [cost, simTraces] = CalyxSim(dM(m,:));
    traces(:,:,m) = simTraces;
    
    if cost < 200
        dataMat(m) = cost;
    else
        dataMat(m) = nan;
    end
    
end

dataMat = [dM, dataMat];
toc
tic
save('CalyxBuildSet.mat','dataMat');
toc
tic
save('CalyxTraces.mat','traces');
toc