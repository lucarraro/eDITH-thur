function [Prod,Conc] = eval_model(betaNew,p0New,tauNew,N_reach,ZCovariateMat,SourceArea,reach_upstream,Qmedian,length_downstream,PathVelocity)
%EVAL_MODEL Summary of this function goes here
%   Detailed explanation goes here
Prod=p0New*exp(ZCovariateMat*betaNew);
    Conc=zeros(N_reach,1);
    for j=1:N_reach
        subset=reach_upstream(reach_upstream(:,j)>0,j);
        Conc(j)=1/Qmedian(j)*sum(Prod(subset).*SourceArea(subset)...
            .*exp(-length_downstream(subset,j)./(tauNew*PathVelocity(subset,j))));
    end
end

