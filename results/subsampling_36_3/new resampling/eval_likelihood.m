function [LoglikNew] = eval_likelihood(rNew,Conc,SitesReach,ReadNo)
%EVAL_POSTERIOR Summary of this function goes here
%   Detailed explanation goes here
    p=1./(Conc(SitesReach)./rNew+1);
    probability=zeros(length(SitesReach),3);
    for i=1:length(SitesReach)
        probability(i,:)=nbinpdf(ReadNo(i,:),rNew,p(i));
    end
    LoglikNew=sum(log(probability(:)));
end

