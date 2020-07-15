function [dists] = determ_vesicle_distances(rand_ves_on_off, num_bins, par)


%rand_ves_on_off variable
    %0: Deterministic 
    %1: Random from distribution
    %5: One SV every 5 nm
    %9: All SVs at 76.8 nm
    %10: All SVs at 95.9 nm (mean of non-integrated Rayleigh)
    %90: SVs at 90 and 110 nm
    %95: SVs at 95 nm
    %122: SVs at 95.9 and 122.1 nm (mean of Rayleigh and integrated Rayleigh)
    %300: SVs at 40, 110, 150, 300 nm
    %1000: 1000 SVs equidist. on [1,400] nm







if sum(rand_ves_on_off == [0 1])
    
    if rand_ves_on_off == 0
        percentiles_step = 1/num_bins;
        percentiles = 0:percentiles_step:1;
        percentiles_dists = (percentiles(1:end-1)+percentiles(2:end))/2;
    %dists = (ray_scale * sqrt(-2*log(1-percentiles_dists)))*1e-3; %calculating quantiles from Rayleigh distr.
    elseif (rand_ves_on_off == 1)
        percentiles_dists = rand(1,num_bins);
    end
    sigm = par(20);
    dists = sqrt( 2 * sigm^2 * gammaincinv(percentiles_dists , 1.5) )*1e-3;


elseif rand_ves_on_off == 9
    dists = ones(1,num_bins)*0.0768;
elseif rand_ves_on_off == 10
    dists = ones(1,num_bins)*0.0958978;   
% elseif rand_ves_on_off == 11
%     dists = ones(1,num_bins)*0.095; 
elseif rand_ves_on_off == 5
    AZ_size = par(16);
    num_dists = abs(AZ_size/(5e-3));
    dists = [1 (5:5:(num_dists*5))]*1e-3;
elseif rand_ves_on_off == 90
    dists = [90 110]*1e-3;
elseif rand_ves_on_off == 95
    dists = ones(1,num_bins)*0.095;
elseif rand_ves_on_off == 122
    dists = [95.9 122.1]*1e-3;    
elseif rand_ves_on_off == 300
    dists = [40 110 150 300]*1e-3;    
elseif rand_ves_on_off == 1000
    dists = (1:399/999:400)*1e-3;
end
