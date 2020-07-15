function [mEPSC_cost] = mEPSC_cost2(pars)

load('mEPSC_aligned.mat')


if (pars(2) <0) || (pars(2) >1)

    mEPSC_cost = 1e6;

else
    

mEPSC_ave_final(mEPSC_ave_final>0) = 0;

time_vec_gen = 0:1e-6:34*1e-3;

mEPSC_gen_func = @(t)(t>=pars(3)).*(pars(1)*(1-exp(-(t-pars(3))/pars(4))).*(pars(2)*exp(-(t-pars(3))/pars(5)) + (1-pars(2)).*exp(-(t-pars(3))/pars(6))));

mEPSC_gen = mEPSC_gen_func(time_vec_gen);

mEPSC_cost_first = ((mEPSC_ave_final - mEPSC_gen).^2)./abs(mEPSC_ave_final);

mEPSC_cost_first(mEPSC_ave_final == 0) = NaN;

mEPSC_cost = sum(mEPSC_cost_first, 'omitnan');

end