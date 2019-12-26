%load ThurData.mat
nrows=size(XX,1); ncols=size(XX,2);

geometry=v2struct(X,Y,XX,YY,Xc,Yc,subcatch,N_reach,AD_pixel,nnodes,outlet,AreaUpstream,reach);

%% Load geology
% load CA (contributing area matrix given outlet position and FD matrix)
CA=imread('DTMad8.tif'); CA=double(CA); CA=flipud(CA); 
CA(CA==-1)=NaN;

% read geology and header lines
fileID = fopen('Clipped_Geology.asc');
GeolOriginal=textscan(fileID,'%f','HeaderLines',6); %f %s %s %f %f %f
fclose(fileID);

% change format and rotate DTM
GeolOriginal=cell2mat(GeolOriginal); GeolOriginal=reshape(GeolOriginal',ncols,nrows);
GeolOriginal=fliplr(GeolOriginal); GeolOriginal=GeolOriginal';

MASK=not(isnan(CA));
Geology=MASK.*GeolOriginal; Geology(Geology==0)=NaN;
clear GeolOriginal

% geology key
% aa=7; ab=4; ac=43; ad=35; ae=97; aj=14; ak=5; al=6; an=2; ao=106; 
% be=3; bf=1; bg=8; bh=9; bi=27; bj=20; hd=89; ie=30; ja=62; jc=57;
% jd=108; je=76; jf=87; jg=48; jh=45; ji=40; jl=37;

% re-attribute classes
Geology(Geology==7)=1001; Geology(Geology==5)=1001; Geology(Geology==2)=1001; % Alluvial 

Geology(Geology==6)=1002; % Moraines 

Geology(Geology==3)=1003; Geology(Geology==1)=1003; Geology(Geology==8)=1003;
Geology(Geology==9)=1003; Geology(Geology==27)=1003; Geology(Geology==20)=1003; % Molasses

Geology(Geology==30)=1004; Geology(Geology==57)=1004; Geology(Geology==108)=1004; 
Geology(Geology==48)=1004; Geology(Geology==45)=1004; % Alpine sediments (penninicum, helveticum)

Geology(Geology==106)=1005; % Water
Geology(Geology==14)=1006; % Loess
Geology(Geology==35)=1007; Geology(Geology==43)=1007; Geology(Geology==97)=1007; % Eboulis/eboulement
Geology(Geology==4)=1008; % Tourbe/peat

%pcolor(Geology); shading flat; colormap('jet'); shading flat; colorbar; title('Geology')

%% load landcover
% read geology and header lines
fileID = fopen('landcover.asc');
LandOriginal=textscan(fileID,'%f','HeaderLines',6); %f %s %s %f %f %f
fclose(fileID); clear text

LandOriginal=cell2mat(LandOriginal); LandOriginal=reshape(LandOriginal',ncols,1801);
LandOriginal=fliplr(LandOriginal); LandOriginal=LandOriginal';
% add upper part to make LandOriginal of the same size of the other matrices
LandOriginal(1802:nrows,:)=NaN;

LandCover=MASK.*LandOriginal; LandCover(LandCover==0)=NaN;
LandCover(LandCover==-9999)=2000;
clear LandOriginal

LandCover(LandCover==1)=2001; % Forest
LandCover(LandCover==2)=2002; LandCover(LandCover==4)=2002; % Rock
LandCover(LandCover==7)=2003; % Urban area
LandCover(LandCover==3)=2004; % Orchard
LandCover(LandCover==9)=2005; % Swamp
LandCover(LandCover==5)=2006; % Lake

%figure; pcolor(XX,YY,LandCover); shading flat; colormap('jet'); shading flat; colorbar; title('Land Cover')

%% Calculate covariates - Land cover
% fraction of local area covered by a given land type
Covariate.LCforest=zeros(N_reach,1); Covariate.LCrock=zeros(N_reach,1);
Covariate.LCurban=zeros(N_reach,1); Covariate.LCorchard=zeros(N_reach,1);
Covariate.LCswamp=zeros(N_reach,1); Covariate.LClake=zeros(N_reach,1);

for i=1:N_reach
    tmp=reach_pixels.(['r',num2str(i)]);
    landcov=LandCover(tmp);
    Covariate.LCforest(i)=sum(landcov==2001)./length(landcov);
    Covariate.LCrock(i)=sum(landcov==2002)./length(landcov);
    Covariate.LCurban(i)=sum(landcov==2003)./length(landcov);
    Covariate.LCorchard(i)=sum(landcov==2004)./length(landcov);
    Covariate.LCswamp(i)=sum(landcov==2005)./length(landcov);
    Covariate.LClake(i)=sum(landcov==2006)./length(landcov);
end

%DrawRiverMap(Covariate.LCrock,1,0,'Fraction of forest','GER',geometry,1,1)

%% Calculate covariates - Geology
% fraction of upstream area covered by a given land type

Covariate.GEOalluvial=zeros(N_reach,1); Covariate.GEOmoraines=zeros(N_reach,1);
Covariate.GEOmolasses=zeros(N_reach,1); Covariate.GEOalpine=zeros(N_reach,1);
Covariate.GEOwater=zeros(N_reach,1); Covariate.GEOloess=zeros(N_reach,1);
Covariate.GEOscree=zeros(N_reach,1); Covariate.GEOpeat=zeros(N_reach,1); 
Covariate.MORupselev=zeros(N_reach,1); Covariate.MORupsslope=zeros(N_reach,1);
for i=1:N_reach
    r_ups=reach_upstream(reach_upstream(:,i)>0,i);
    tmp=[];
    for ind_r=1:length(r_ups)
        tmp=[tmp; reach_pixels.(['r',num2str(r_ups(ind_r))])];
    end
    geol=Geology(tmp); 
    Covariate.GEOalluvial(i)=sum(geol==1001)./length(geol);
    Covariate.GEOmoraines(i)=sum(geol==1002)./length(geol);
    Covariate.GEOmolasses(i)=sum(geol==1003)./length(geol);
    Covariate.GEOalpine(i)=sum(geol==1004)./length(geol);
    Covariate.GEOwater(i)=sum(geol==1005)./length(geol);
    Covariate.GEOloess(i)=sum(geol==1006)./length(geol);
    Covariate.GEOscree(i)=sum(geol==1007)./length(geol);
    Covariate.GEOpeat(i)=sum(geol==1008)./length(geol);
    Covariate.MORupselev(i)=sum(LocalElevation(r_ups).*AreaLocal(r_ups)/AreaUpstream(i));
    Covariate.MORupsslope(i)=sum(reach_slope(r_ups).*length_reach(r_ups)/sum(length_reach(r_ups)));
end

Covariate.MORarea=AreaUpstream;
Covariate.MORlocalslope=reach_slope;
Covariate.MORlocalelev=LocalElevation;
Covariate.MORstreamorder=stream_order_reach;
%Covariate.MORdisttooutlet=length_downstream(:,find(down_reach==0)); 
%If I add DistanceToOutlet as covariate, VIF of local elev goes >10 

%% Evaluate multicollinearity
CovariateReduced=Covariate;
% remove covariates with highest VIF until all remaining covariates have VIF<10
CovariateReduced=rmfield(CovariateReduced,{'GEOmolasses','MORupselev'}); 
CovariateTable=struct2table(CovariateReduced); 
CovariateMat=cell2mat(struct2cell(CovariateReduced)');

R=corrcoef(CovariateMat);
VIF=diag(inv(R))'

% z-transform
ZCovariateMat=zeros(size(CovariateMat));
for i=1:size(CovariateMat,2)
    ZCovariateMat(:,i)=(CovariateMat(:,i)-mean(CovariateMat(:,i)))./std(CovariateMat(:,i));
end
CovariateNames=fieldnames(CovariateReduced);
%DrawRiverMap(ZCovariateMat(:,6),max(ZCovariateMat(:,6)),min(ZCovariateMat(:,6)),'Fraction of loess','GER',geometry,1,1)

save('data/Covariates.mat','ZCovariateMat','CovariateMat','CovariateNames')