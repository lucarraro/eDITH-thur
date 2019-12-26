
SumGenusPerSite=sum(GenusPerSite,2);

RegioCov;
ZCovMat=[ZCovariateMat RegioCovMat(:,1:end-1)];
CovNames=[CovariateNames; RegioCovNames];

%% calculate lateral Q and source area
Qlat=zeros(N_reach,1);
for i=1:N_reach
    subset=find(down_reach==i);
    if ~isempty(subset)
        Qlat(i)=Qmedian(i)-sum(Qmedian(subset));
    else
        Qlat(i)=Qmedian(i);
    end
end
SourceArea=ReachWidth.*length_reach;
PathVelocity=zeros(N_reach);
for i=1:N_reach
    for j=1:N_reach
        path=list_reach_downstream{i,j};
        if ~isempty(path)
            PathVelocity(i,j)=length_downstream(i,j)/(sum(length_reach(path)./VelocityMedian(path)));
        end
    end
end

N_param=size(ZCovMat,2)+2; 


DecayTime=zeros(length(GenusName),1); CovariateSign=zeros(N_param-2,length(GenusName));
PresenceMat=zeros(N_reach,length(GenusName)); DetectionProbAll=zeros(N_reach,length(GenusName));

for g=1:length(GenusName)
    Genus=GenusName{g};  
    load(['../results_MCMC_all_NewCovSet/',Genus])
    quant_par=zeros(N_param,3);
    for i=1:N_param
        quant_par(i,:)=quantile(par(:,i),[0.025 0.5 0.975]);
    end
    ResultsAll.(Genus).DecayTime=exp(quant_par(end,2))/3600;
    ResultsAll.(Genus).param=quant_par(:,2);
    ResultsAll.(Genus).prod=quant_p(12,:)';
    ResultsAll.(Genus).Conc=quant_C(12,:)';
    
    UnconnectedConc=quant_p(12,:)'.*length_reach.*ReachWidth./Qlat.*exp(-length_reach./VelocityMedian./(exp(quant_par(end,2))));
    ResultsAll.(Genus).DetectionProbAll=UnconnectedConc./(1+UnconnectedConc);
    ResultsAll.(Genus).Presence= ResultsAll.(Genus).DetectionProbAll > 2/3;
    ResultsAll.(Genus).CovariateSign=(quant_par(1:N_param-2,1)>0) - (quant_par(1:N_param-2,3)<0);
    
    DecayTime(g)=ResultsAll.(Genus).DecayTime;
    CovariateSign(:,g)=ResultsAll.(Genus).CovariateSign;
    DetectionProbAll(:,g)=ResultsAll.(Genus).DetectionProbAll; 
    PresenceMat(:,g)=ResultsAll.(Genus).Presence;  
    
    for indSite=1:numel(SitesReach)
        ResultsAll.(Genus).p_value(indSite,1)=GOF_NBtest(ReadNumbers.(Genus)(indSite,:)',ResultsAll.(Genus).Conc(SitesReach(indSite)),100000);
    end
    AcceptedAllSites(g,1)=sum(ResultsAll.(Genus).p_value>0.05);
    AcceptedPerSite(g,:)=ResultsAll.(Genus).p_value>0.05;
end

%% Kicknet vs eDNA
KicknetSiteReach=zeros(length(KicknetSiteID),1);
for i=1:length(KicknetSiteID)
    k=KicknetSiteID(i);
    KicknetSiteReach(i)=SitesReach(find(siteID==k));
end
KicknetPresence=KicknetData>0;

Site_eDNAkick=find(ismember(siteID,KicknetSiteID));% ReadNumbers.(Genus)(Site_eDNAKick) are read numbers at the same sites as kicknet

KicknetVsModel_mat=zeros(length(KicknetSiteID),length(GeneraOrder));
Scores_eDNAKick=zeros(length(GeneraOrder),1);
TruePositiveAll=zeros(length(GeneraOrder),1); TrueNegativeAll=zeros(length(GeneraOrder),1);
FalsePositiveAll=zeros(length(GeneraOrder),1); FalseNegativeAll=zeros(length(GeneraOrder),1);
TruePositive_eDNA=zeros(length(GeneraOrder),1); TrueNegative_eDNA=zeros(length(GeneraOrder),1);
FalsePositive_eDNA=zeros(length(GeneraOrder),1); FalseNegative_eDNA=zeros(length(GeneraOrder),1);

for i=1:length(GeneraOrder)
    Genus=GeneraOrder{i,1};
    % read model
    if sum(strcmp(cellstr(GenusName),Genus))>0
        ind_model=find(strcmp(cellstr(GenusName),Genus));
        presence_probAll=DetectionProbAll(KicknetSiteReach,ind_model);
        presence_probAll(presence_probAll<1e-50)=1e-50;

        % read eDNA data
        presence_eDNA = sum(ReadNumbers.(Genus)>0,2);    
    else
        presence_probAll=1e-50*ones(length(KicknetSiteReach),1);
        presence_eDNA = zeros(length(siteID),1);
    end
    % read kicknet data
    if sum(strcmp(cellstr(KicknetName),Genus))>0
        ind_kicknet=find(strcmp(cellstr(KicknetName),Genus));
        presence_kick=KicknetPresence(:,ind_kicknet);
    else
        presence_kick=zeros(length(KicknetSiteReach),1);
    end
    % calculate probability
    for j=1:length(KicknetSiteReach)
        if presence_kick(j)==1
            Scores_eDNAKick(i) = Scores_eDNAKick(i) + log(presence_probAll(j));
            % model
            if presence_probAll(j) > 2/3
                TruePositiveAll(i) = TruePositiveAll(i) + 1;
                KicknetVsModel_mat(j,i)=1; 
            else
                FalseNegativeAll(i) = FalseNegativeAll(i) + 1;
                KicknetVsModel_mat(j,i)=2; 
            end
            % eDNA measurements
            if presence_eDNA(Site_eDNAkick(j)) > 1
                TruePositive_eDNA(i) = TruePositive_eDNA(i) + 1;
            else
                FalseNegative_eDNA(i) = FalseNegative_eDNA(i) + 1;
            end
            
        else
            Scores_eDNAKick(i) = Scores_eDNAKick(i) + log(1-presence_probAll(j));
            if presence_probAll(j) > 2/3
                FalsePositiveAll(i) = FalsePositiveAll(i) + 1;
                KicknetVsModel_mat(j,i)=3; 
            else
                TrueNegativeAll(i) = TrueNegativeAll(i) + 1;
                KicknetVsModel_mat(j,i)=4; 
            end
            % eDNA measurements
            if presence_eDNA(Site_eDNAkick(j)) > 1
                FalsePositive_eDNA(i) = FalsePositive_eDNA(i) + 1;
            else
                TrueNegative_eDNA(i) = TrueNegative_eDNA(i) + 1;
            end
        end
    end
end
AccuracyAll=(TruePositiveAll+TrueNegativeAll)/length(KicknetSiteReach);
Accuracy_eDNA=(TruePositive_eDNA+TrueNegative_eDNA)/length(KicknetSiteReach);

% display accuracy by TP, TN, FP, FN
[~,indFN]=sort(1e6*FalseNegativeAll-TruePositiveAll-TrueNegativeAll,'ascend');
indFN=indFN(ismember(GeneraOrder(indFN,1),GenusName)); % pick only the genera detected by eDNA