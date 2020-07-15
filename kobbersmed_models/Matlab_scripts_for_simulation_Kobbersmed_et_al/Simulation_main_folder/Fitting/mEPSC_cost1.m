function [mEPSC_cost] = mEPSC_cost1(pars)
% 
% load('mEPSC_new.mat')
% 
% time_vec_gen = 0:1e-6:34*1e-3;
% 
% mEPSC_gen_func = @(t)(t>=pars(3)).*(pars(1)*(1-exp(-(t-pars(3))/pars(4))).*(pars(2)*exp(-(t-pars(3))/pars(5)) + (1-pars(2)).*exp(-(t-pars(3))/pars(6))));
% 
% mEPSC_gen = mEPSC_gen_func(time_vec_gen);
% 
% mEPSC_cost_first = ((mEPSC_current(1:end-1) - mEPSC_gen(1:end-1)).^2)./abs(mEPSC_current(1:end-1));
% 
% % mEPSC_cost_first(mEPSC_cost_first == inf) = NaN;
% 
% mEPSC_cost = sum(mEPSC_cost_first, 'omitnan');
% 


% disp('Before load')



if (pars(2) <0) || (pars(2) >1)

    mEPSC_cost = 1e6;

else
    
    load('mEPSC_new.mat')
    
%     disp('After load')

%     mEPSC_ave(mEPSC_ave>0) = 0;
    
    
    peak_current = min(mEPSC_current);
    ind_10perc = find(mEPSC_current<0.1*peak_current, 1, 'first');
    ind_34ms = find(mEPSC_time==0.034);

    time_vec_gen = 0:1e-6:34*1e-3;

    mEPSC_gen_func = @(t)(t>=pars(3)).*(pars(1)*(1-exp(-(t-pars(3))/pars(4))).*(pars(2)*exp(-(t-pars(3))/pars(5)) + (1-pars(2)).*exp(-(t-pars(3))/pars(6))));

    mEPSC_gen = mEPSC_gen_func(time_vec_gen);

    mEPSC_cost_first = ((mEPSC_current(ind_10perc:ind_34ms) - mEPSC_gen(ind_10perc:ind_34ms)).^2)./abs(mEPSC_current(ind_10perc:ind_34ms));

%     mEPSC_cost_first(mEPSC_ave_final == 0) = NaN;

    mEPSC_cost = sum(mEPSC_cost_first, 'omitnan');
    
%     disp('After cost calculation')

end