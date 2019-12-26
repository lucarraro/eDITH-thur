function p = GOF_NBtest(data,Conc,N_sample)
% GOF test for NB
% according to Mi et al., Goodness-of-Fit tests and model diagnostics for
% negative binomial regression of RNA sequencing data, PLOS one, 2015
p=1/(1+Conc);

mean_NB=Conc;
s_NB=1e-10+sqrt(sum((data-mean_NB).^2));
res_data=(data-mean_NB)/s_NB;

MCdata=nbinrnd(1,p,3,N_sample);
s_MC=1e-10+sqrt(sum((MCdata-mean_NB).^2));
res_MCdata=(MCdata-mean_NB)./s_MC;

d_0=sum((res_data-median(res_MCdata')').^2);
d_h=sum((res_MCdata-median(res_MCdata')').^2);

p=(sum(d_h>=d_0)+1)/(N_sample+1);