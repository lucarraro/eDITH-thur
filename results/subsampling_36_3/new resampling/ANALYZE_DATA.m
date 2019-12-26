clear all; close all; clc

load GenusData
load Covariates
load ThurData
load ThurHydrology

CovariateNames={'L-FO';'L-RO';'L-UR';'L-OR';'L-SW';'L-LA';'G-AL';'G-MO';'G-AP';'G-WA';
    'G-LO';'G-SC';'G-PE';'M-US';'M-DA';'M-LS';'M-LE';'M-SO'};

%% Table with ID and genera names

for i=1:length(GeneraOrder)
    tmp{i,1}=[GeneraOrder{i,2},GeneraOrder{i,1}];
end
[~,indAlphabetAll]=sort(tmp);

% find genera keys
for i=1:length(GenusName)
    eDNA_to_all(i,1)=find(ismember(GeneraOrder(:,1),GenusName(i)));
end
for i=1:length(GeneraOrder)
    all_to_alphabet(i,1)=find(indAlphabetAll==i);
end
all_to_eDNA=nan(length(GeneraOrder),1);
for i=1:length(GeneraOrder)
    tmp=find(ismember(GenusName,GeneraOrder(i,1)));
    if not(isempty(tmp))
        all_to_eDNA(i,1)=tmp;
    end
end

try load('results_All.mat')
catch
    read_ALL_results;
    save('results_All.mat','AcceptedAllSites','AcceptedPerSite','AccuracyAll','indFN','N_param','Qlat')
end

%% READ Results
try load('results_valid_NewResampling.mat')
catch
    Results=[];
    for N_ValidSites=12:12:36
        for RUN=1:3
            N_ValidSites
            RUN
            Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).AcceptedAllSites=zeros(length(GenusName),1);
            Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).DetectionProb=zeros(length(GenusName),N_reach);
            for g=1:length(GenusName)
                Genus=GenusName{g};
                load(['results_valid_',num2str(N_ValidSites),'_',num2str(RUN),'/',Genus])
                quant_par=zeros(N_param,3);
                for i=1:N_param
                    quant_par(i,:)=quantile(par(:,i),[0.025 0.5 0.975]);
                end
                UnconnectedConc=quant_p(12,:)'.*length_reach.*ReachWidth./Qlat.*exp(-length_reach./VelocityMedian./(exp(quant_par(end,2))));
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).DetectionProb=UnconnectedConc./(1+UnconnectedConc);
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).ValidSites=ValidSites;
                CalibSites=setdiff(1:61,ValidSites);
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).Conc=quant_C(12,:)';
                for indSite=1:numel(SitesReach)
                    Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).p_value(indSite,1)=...
                        GOF_NBtest(ReadNumbers.(Genus)(indSite,:)',Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).Conc(SitesReach(indSite)),100000);
                end
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).GOF(g,1)=...
                    sum(Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).p_value>0.05);
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).GOFcalib(g,1)=...
                    sum(Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).p_value(CalibSites)>0.05);
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).GOFvalid(g,1)=...
                    sum(Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).p_value(ValidSites)>0.05);
                Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).DetectionProb(g,:)=...
                    Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).(Genus).DetectionProb;
            end
            [TP,TN,FP,FN] = eval_accuracy(Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).DetectionProb,...
                KicknetPresence,Site_eDNAkick,SitesReach,2/3,GenusName,KicknetName);
            Results.(['v',num2str(N_ValidSites)]).(['r',num2str(RUN)]).Accuracy=(TP+TN)/length(Site_eDNAkick);
        end
    end
    ValidSites_v12=[Results.v12.r1.ValidSites Results.v12.r2.ValidSites Results.v12.r3.ValidSites];
    ValidSites_v24=[Results.v24.r1.ValidSites Results.v24.r2.ValidSites Results.v24.r3.ValidSites];
    ValidSites_v36=[Results.v36.r1.ValidSites Results.v36.r2.ValidSites Results.v36.r3.ValidSites];   
    GOF_v12=[Results.v12.r1.GOF(all_to_eDNA(indFN)) Results.v12.r2.GOF(all_to_eDNA(indFN)) Results.v12.r3.GOF(all_to_eDNA(indFN))];
    GOF_v24=[Results.v24.r1.GOF(all_to_eDNA(indFN)) Results.v24.r2.GOF(all_to_eDNA(indFN)) Results.v24.r3.GOF(all_to_eDNA(indFN))];
    GOF_v36=[Results.v36.r1.GOF(all_to_eDNA(indFN)) Results.v36.r2.GOF(all_to_eDNA(indFN)) Results.v36.r3.GOF(all_to_eDNA(indFN))];
    GOFcalib_v12=[Results.v12.r1.GOFcalib(all_to_eDNA(indFN)) Results.v12.r2.GOFcalib(all_to_eDNA(indFN)) Results.v12.r3.GOFcalib(all_to_eDNA(indFN))];
    GOFcalib_v24=[Results.v24.r1.GOFcalib(all_to_eDNA(indFN)) Results.v24.r2.GOFcalib(all_to_eDNA(indFN)) Results.v24.r3.GOFcalib(all_to_eDNA(indFN))];
    GOFcalib_v36=[Results.v36.r1.GOFcalib(all_to_eDNA(indFN)) Results.v36.r2.GOFcalib(all_to_eDNA(indFN)) Results.v36.r3.GOFcalib(all_to_eDNA(indFN))];
    GOFvalid_v12=[Results.v12.r1.GOFvalid(all_to_eDNA(indFN)) Results.v12.r2.GOFvalid(all_to_eDNA(indFN)) Results.v12.r3.GOFvalid(all_to_eDNA(indFN))];
    GOFvalid_v24=[Results.v24.r1.GOFvalid(all_to_eDNA(indFN)) Results.v24.r2.GOFvalid(all_to_eDNA(indFN)) Results.v24.r3.GOFvalid(all_to_eDNA(indFN))];
    GOFvalid_v36=[Results.v36.r1.GOFvalid(all_to_eDNA(indFN)) Results.v36.r2.GOFvalid(all_to_eDNA(indFN)) Results.v36.r3.GOFvalid(all_to_eDNA(indFN))];
    Accuracy_v12=[Results.v12.r1.Accuracy(all_to_eDNA(indFN)) Results.v12.r2.Accuracy(all_to_eDNA(indFN)) Results.v12.r3.Accuracy(all_to_eDNA(indFN))];
    Accuracy_v24=[Results.v24.r1.Accuracy(all_to_eDNA(indFN)) Results.v24.r2.Accuracy(all_to_eDNA(indFN)) Results.v24.r3.Accuracy(all_to_eDNA(indFN))];
    Accuracy_v36=[Results.v36.r1.Accuracy(all_to_eDNA(indFN)) Results.v36.r2.Accuracy(all_to_eDNA(indFN)) Results.v36.r3.Accuracy(all_to_eDNA(indFN))];
    save('results_valid_NewResampling.mat','GOF_v12','GOF_v24','GOF_v36','GOFcalib_v12','GOFcalib_v24','GOFcalib_v36',...
        'GOFvalid_v12','GOFvalid_v24','GOFvalid_v36','Accuracy_v12','Accuracy_v24','Accuracy_v36',...
        'ValidSites_v12','ValidSites_v24','ValidSites_v36')
end

%% figure GOF
figure('units','centimeters','position',[0 0 25 10]); 
bar(AcceptedAllSites(all_to_eDNA(indFN))/61,'FaceColor','n'); hold on; box off
set(gca,'tickdir','out','ytick',[0:0.2:1],'xtick',[1:50],'ylim',[0.6 1])
ylabel('Fraction of sites where H_0 cannot be rejected')

plot((1:50)-0.2,GOF_v12(:,1)/61,'ob');
plot((1:50)-0.2,GOF_v12(:,2)/61,'+b');
plot((1:50)-0.2,GOF_v12(:,3)/61,'xb');
for i=1:length(GenusName)
    plot([i-0.2 i-0.2],[min(GOF_v12(i,:)) max(GOF_v12(i,:))]/61,'b')
end

plot((1:50),GOF_v24(:,1)/61,'or');
plot((1:50),GOF_v24(:,2)/61,'+r');
plot((1:50),GOF_v24(:,3)/61,'xr');
for i=1:length(GenusName)
    plot([i i],[min(GOF_v24(i,:)) max(GOF_v24(i,:))]/61,'r')
end

plot((1:50)+0.2,GOF_v36(:,1)/61,'og');
plot((1:50)+0.2,GOF_v36(:,2)/61,'+g');
plot((1:50)+0.2,GOF_v36(:,3)/61,'xg');
for i=1:length(GenusName)
    plot([i i]+0.2,[min(GOF_v36(i,:)) max(GOF_v36(i,:))]/61,'g')
end


%% note that LossOfGOF is already well ordered
LossOfGOF=zeros(length(GenusName),3);
for g=1:length(GenusName)
    LossOfGOF(g,1)=mean(GOF_v12(g,:)-AcceptedAllSites(all_to_eDNA(indFN(g))))/61;
    LossOfGOF(g,2)=mean(GOF_v24(g,:)-AcceptedAllSites(all_to_eDNA(indFN(g))))/61;
    LossOfGOF(g,3)=mean(GOF_v36(g,:)-AcceptedAllSites(all_to_eDNA(indFN(g))))/61;
end

LossOfGOFcalib=zeros(length(GenusName),3);
for g=1:length(GenusName)
    tmp12=zeros(1,3); tmp24=zeros(1,3); tmp36=zeros(1,3);
    for r=1:3
      tmp12(:,r)=sum(AcceptedPerSite(all_to_eDNA(indFN(g)),setdiff(1:61,ValidSites_v12(:,r))));
      tmp24(:,r)=sum(AcceptedPerSite(all_to_eDNA(indFN(g)),setdiff(1:61,ValidSites_v24(:,r))));
      tmp36(:,r)=sum(AcceptedPerSite(all_to_eDNA(indFN(g)),setdiff(1:61,ValidSites_v36(:,r))));
    end
    LossOfGOFcalib(g,1)=mean(GOFcalib_v12(g,:)-tmp12)/(61-12);
    LossOfGOFcalib(g,2)=mean(GOFcalib_v24(g,:)-tmp24)/(61-24);
    LossOfGOFcalib(g,3)=mean(GOFcalib_v36(g,:)-tmp36)/(61-36);
end

LossOfGOFvalid=zeros(length(GenusName),3);
for g=1:length(GenusName)
    tmp12=zeros(1,3); tmp24=zeros(1,3); tmp36=zeros(1,3);
    for r=1:3
      tmp12(:,r)=sum(AcceptedPerSite(all_to_eDNA(indFN(g)),ValidSites_v12(:,r)));
      tmp24(:,r)=sum(AcceptedPerSite(all_to_eDNA(indFN(g)),ValidSites_v24(:,r)));
      tmp36(:,r)=sum(AcceptedPerSite(all_to_eDNA(indFN(g)),ValidSites_v36(:,r)));
    end
    LossOfGOFvalid(g,1)=mean(GOFvalid_v12(g,:)-tmp12)/12;
    LossOfGOFvalid(g,2)=mean(GOFvalid_v24(g,:)-tmp24)/24;
    LossOfGOFvalid(g,3)=mean(GOFvalid_v36(g,:)-tmp36)/36;
end

% mean (across genera) Loss of accuracy due to reduction in number of sampling sites
MeanLossOfGOF = mean(LossOfGOF)
MeanLossOfGOFcalib = mean(LossOfGOFcalib)
MeanLossOfGOFvalid = mean(LossOfGOFvalid)

% ANOVA
anova1(([GOF_v12 GOF_v24 GOF_v36]-AcceptedAllSites(all_to_eDNA(indFN)))/61,{'v12','v12','v12','v24','v24','v24','v36','v36','v36'})

%% figure accuracy
figure('units','centimeters','position',[0 0 25 10]); 
bar(AccuracyAll(indFN),'FaceColor','n'); hold on; box off
set(gca,'tickdir','out','ytick',[0.25:0.25:1],'xtick',[1:50],'ylim',[0.25 1])
ylabel('Accuracy')

plot((1:50)-0.2,Accuracy_v12(:,1),'ob');
plot((1:50)-0.2,Accuracy_v12(:,2),'+b');
plot((1:50)-0.2,Accuracy_v12(:,3),'xb');
for i=1:length(GenusName)
    plot([i-0.2 i-0.2],[min(Accuracy_v12(i,:)) max(Accuracy_v12(i,:))],'b')
end

plot((1:50),Accuracy_v24(:,1),'or');
plot((1:50),Accuracy_v24(:,2),'+r');
plot((1:50),Accuracy_v24(:,3),'xr');
for i=1:length(GenusName)
    plot([i i],[min(Accuracy_v24(i,:)) max(Accuracy_v24(i,:))],'r')
end

plot((1:50)+0.2,Accuracy_v36(:,1),'og');
plot((1:50)+0.2,Accuracy_v36(:,2),'+g');
plot((1:50)+0.2,Accuracy_v36(:,3),'xg');
for i=1:length(GenusName)
    plot([i i]+0.2,[min(Accuracy_v36(i,:)) max(Accuracy_v36(i,:))],'g')
end

%% note that LossOfAccuracy is already well ordered
LossOfAccuracy=zeros(length(GenusName),3);
for g=1:length(GenusName)
    LossOfAccuracy(g,1)=mean(Accuracy_v12(g,:)-AccuracyAll(indFN(g)));
    LossOfAccuracy(g,2)=mean(Accuracy_v24(g,:)-AccuracyAll(indFN(g)));
    LossOfAccuracy(g,3)=mean(Accuracy_v36(g,:)-AccuracyAll(indFN(g)));
end
% mean (across genera) Loss of accuracy due to reduction in number of sampling sites
MeanLossOfAccuracy = mean(LossOfAccuracy)

% ANOVA
anova1(([Accuracy_v12 Accuracy_v24 Accuracy_v36]-AccuracyAll(indFN)),{'v12','v12','v12','v24','v24','v24','v36','v36','v36'})
