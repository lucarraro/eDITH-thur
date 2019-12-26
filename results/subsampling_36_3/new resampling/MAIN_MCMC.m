clear all; close all; clc

rng('shuffle')

CORE = 1; N_CORES = 1; RUN = 3; N_ValidSites = 36;

load ThurData
load ThurHydrology
load Covariates
load GenusData

geometry=v2struct(X,Y,XX,YY,Xc,Yc,subcatch,N_reach,AD_pixel,nnodes,outlet,AreaUpstream,reach);

PathVelocity=zeros(N_reach);
for i=1:N_reach
    for j=1:N_reach
        path=list_reach_downstream{i,j};
        if ~isempty(path)
            PathVelocity(i,j)=length_downstream(i,j)/(sum(length_reach(path)./VelocityMedian(path)));
        end
    end
end
SourceArea=ReachWidth.*length_reach;
quant_vec=[0:0.005:0.05 0.075:0.025:0.125 0.15:0.05:0.85 0.875:0.025:0.925 0.95:0.005:1];

%% define validation sites
RareSpecies=find(sum(GenusPerSite,2)<=10);
NeverUndetected=sum(GenusPerSite(RareSpecies,:)); NeverUndetected=find(NeverUndetected==0); % sites where rare species were never detected
%stream_order_reach(SitesReach(tmp))
for i=1:length(GenusName)
    Genus=GenusName{i};
    [a,b]=find(ReadNumbers.(Genus)==max(max(ReadNumbers.(Genus))));
    SiteMaxRead(i,1)=a;
end
NeverMax=setdiff(1:61,SiteMaxRead);
intersect(NeverUndetected,NeverMax); % feasible sites

% find how many valid sites per each stream order
N_CalibSites=length(siteX)-N_ValidSites;
tableDHondt=zeros(length(siteX),5);
for i=1:5
    n_ReachPerStreamOrder(1,i)=sum(stream_order_reach(SitesReach)==i);
    for j=1:length(siteX)
        table_DHondt(j,i)=n_ReachPerStreamOrder(1,i)/j;
    end
end
DistribAllSites=histcounts(stream_order_reach(SitesReach));
vector_DHondt=table_DHondt(:);
[vector_DHondt,indices_DHondt]=sort(vector_DHondt,'descend');
indices_DHondt=floor((indices_DHondt-1)/length(siteX))+1;
DistribCalibSites=histcounts(indices_DHondt(1:N_CalibSites));
DistribValidSites=DistribAllSites-DistribCalibSites;

% indicize SiteReach by stream order
SortbyStreamOrder=zeros(length(siteX),1);
tmp=zeros(5,1);
for i=1:length(siteX)    
    k=tmp(stream_order_reach(SitesReach(i)));
    SortbyStreamOrder(i)=k+1;
    tmp(stream_order_reach(SitesReach(i)))=k+1;
end

% sample validation sites according to the previous partitioning
ValidSites=[];
for i=1:5
    tmp=randsample(find(stream_order_reach(SitesReach)==i),DistribValidSites(i));
    ValidSites=[ValidSites; tmp];
end

ValidSites % pick 12 among the above feasible sites (check for homogeneity in geographical distribution and stream order values) % 2 4 10 18 20 24 28 31 36 40 42 60
CalibSites=setdiff(1:length(SitesReach),ValidSites);

%% Choose covariates
choose_covariates=[1:35];
burnin=1000; chain_length=5000;

% choose type of covariates
%ZCovMat=ZCovariateMat(:,choose_covariates);
RegioCov;
ZCovMat=[ZCovariateMat RegioCovMat(:,1:end-1)];

%CovNames=CovariateNames(choose_covariates);

N_run=5e7; N_param=length(choose_covariates)+2;

q = @(COVmat,N_param) mvnrnd(zeros(1,N_param),COVmat);

vec_indGenus = randperm(length(GenusName));
for indGenus_noShuffle = 1:length(GenusName)
    indGenus=vec_indGenus(indGenus_noShuffle);
    tic
    COVmatReal=0.1*eye(N_param);
    Genus=GenusName{indGenus};   
    Loglik=nan(burnin+chain_length,1); par=zeros(burnin+chain_length,N_param);
    Prod_mat=zeros(burnin+chain_length,N_reach); Conc_mat=zeros(burnin+chain_length,N_reach); 
    % find initial parameter set
    LoglikInit=-inf(500,1); parInit=zeros(N_param,500);
    disp(sprintf('%s: Searcing for initial parameter set...',Genus))
    k=0;
    while max(LoglikInit)==-Inf
        k=k+1;
        disp(sprintf('round %d',k))
    for indPart = 1:500
        parInit(:,indPart)=[normrnd(0,3,N_param-2,1); normrnd(-10,3); normrnd(9,0.5)]; % initialize parameters' values
        [Prod,Conc] = eval_model(parInit(1:N_param-2,indPart),exp(parInit(N_param-1,indPart)),exp(parInit(N_param,indPart)),...
            N_reach,ZCovMat,SourceArea,reach_upstream,Qmedian,length_downstream,PathVelocity);
        LoglikInit(indPart,1) = sum(log(normpdf(parInit([1:N_param-2 N_param],indPart),[zeros(N_param-2,1); 9],[3*ones(N_param-2,1); 0.5]))) + ...
            eval_likelihood(1,Conc,SitesReach(CalibSites),ReadNumbers.(Genus)(CalibSites,:));
        
    end
    end   
    indAccept=1; indInit=find(LoglikInit==max(LoglikInit));
    par(1,:)=parInit(:,indInit); Loglik(1)=LoglikInit(indInit);
    [Prod,Conc] = eval_model(par(1:N_param-2)',exp(par(N_param-1)),exp(par(N_param)),...
            N_reach,ZCovMat,SourceArea,reach_upstream,Qmedian,length_downstream,PathVelocity);
    prod_mat(1,:)=Prod'; Conc_mat(1,:)=Conc'; 
    parOld=par(1,:); LoglikOld=Loglik(1);
    disp(sprintf('%s: Accepted: %d  -  Total: 1  -  Loglik %.1f  - tau: %.1f h  -  Elapsed time: %.0f s',...
                GenusName{indGenus},indAccept,LoglikOld,exp(parOld(N_param))/3600,toc));
    rej=0;
    for indRun = 2:N_run
        if mod(indAccept,10)==0
            COVmatReal=((2.38/sqrt(N_param))^2 * cov(par(1:indAccept-1,1:N_param)) + 1e-4 * eye(N_param));
        end
        COVmat=COVmatReal*exp(-rej/5000);  
        parNew = parOld + q(COVmat,N_param);
        [Prod,Conc] = eval_model(parNew(1:N_param-2)',exp(parNew(N_param-1)),exp(parNew(N_param)),...
            N_reach,ZCovMat,SourceArea,reach_upstream,Qmedian,length_downstream,PathVelocity);
        LoglikNew = sum(log(normpdf(parNew([1:N_param-2 N_param]),[zeros(1,N_param-2) 9],[3*ones(1,N_param-2) 0.5]))) + ...
            eval_likelihood(1,Conc,SitesReach(CalibSites),ReadNumbers.(Genus)(CalibSites,:));
        if LoglikNew > LoglikOld || rand<exp(LoglikNew-LoglikOld)
            rej=0;
            indAccept=indAccept+1; 
            parOld=parNew; par(indAccept,:)=parNew;
            LoglikOld=LoglikNew; Loglik(indAccept)=LoglikNew;
            Prod_mat(indAccept,:)=Prod'; Conc_mat(indAccept,:)=Conc';
            disp(sprintf('No. ValidSites: %d  -  RUN: %d  -  %s: Accepted: %d  -  Total: %d  -  Loglik %.1f  - tau: %.1f h  -  Elapsed time: %.0f s',...
                N_ValidSites,RUN,GenusName{indGenus},indAccept,indRun,LoglikOld,exp(parOld(N_param))/3600,toc));
        else rej=rej+1; % reduce covariance matrix if the previous set is rejected             
        end            
        if indAccept==burnin+chain_length
            break
        end
    end
    % cut burn-in
    Loglik(1:burnin)=[]; par(1:burnin,:)=[]; Prod_mat(1:burnin,:)=[]; Conc_mat(1:burnin,:)=[];
    % Prod & Conc quantiles
    vec_quantiles=[0 0.005 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85 0.9 0.95 0.975 0.99 0.995 1];
    quant_p=quantile(Prod_mat,vec_quantiles); 
    quant_C=quantile(Conc_mat,vec_quantiles); 
    save(['results_valid_',num2str(N_ValidSites),'_',num2str(RUN),'/',Genus],'par','Loglik','quant_p','quant_C','vec_quantiles','ValidSites')
    disp(' ')
end
% R2_avg=1-sum((Conc(SitesReach)'-mean(ReadNumbers.(Genus)')).^2)/sum((mean(ReadNumbers.(Genus)')-mean(mean(ReadNumbers.(Genus)'))).^2);


