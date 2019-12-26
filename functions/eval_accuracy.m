function [TP,TN,FP,FN] = eval_accuracy(DetectionProb,KicknetPresence,Site_eDNAkick,SitesReach,threshold,GenusName,KicknetName)

% DetectionProb is matrix for a given model realization (all genera x all reaches)
% KicknetPresence is matrix (all genera x sampled reaches)
% (GeneraOrder(find(ismember(GeneraOrder(:,1),GenusName(k,1))),1) = GenusName(k,1)) 
% if genus is absent from kicknet, then KicknetPresence = vector of zeros

% site i (kicknet) -> site SitesReach(Site_eDNAkick(i)) (eDNA)

TP=zeros(length(GenusName),1); TN = zeros(length(GenusName),1); FP = zeros(length(GenusName),1); FN = zeros(length(GenusName),1);

for g=1:length(GenusName)
    % find vector of kicknet presence for given genus
    ind = find(ismember(KicknetName(:,1),GenusName(g,1)));
    if isempty(ind)
        Kicknet_pres=zeros(1,length(Site_eDNAkick));
    else
        Kicknet_pres=KicknetPresence(:,ind)';
    end
    
    for i=1:length(Site_eDNAkick)
        j=SitesReach(Site_eDNAkick(i));
        if DetectionProb(g,j) >= threshold && Kicknet_pres(i) == 1
            TP(g) = TP(g) + 1;
        elseif DetectionProb(g,j) < threshold && Kicknet_pres(i) == 0
            TN(g) = TN(g) + 1;
        elseif DetectionProb(g,j) >= threshold && Kicknet_pres(i) == 0
            FP(g) = FP(g) + 1;
        elseif DetectionProb(g,j) < threshold && Kicknet_pres(i) == 1
            FN(g) = FN(g) + 1;
        end
    end
    
end
    
    
