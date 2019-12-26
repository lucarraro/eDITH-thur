RegioCovNames={'TH1'; %Thur I';
               'TH2'; %Thur II';
               'TH3'; %'Thur III';
               'TH4'; %'Thur IV';
               'TH5'; %'Thur V';
               'TH6'; %'Thur VI';
               'TH7'; %'Thur VII';
               'TH8'; %'Thur VIII';
               'GL1'; %Glatt I';
               'GL2'; %'Glatt II';
               'WIS'; %Wissbach';
               'NE1'; %Necker I';
               'NE2'; %'Necker II';
               'NE3'; %'Necker III';
               'GON'; %Gonzenbach';
               'DIE';%Dietfurterbach';
               'LUT'}; %'Luteren'
               
           
RegioCovID=ones(760,1); % Thur I
RegioCovID(reach_upstream(reach_upstream(:,312)>0,312))=2; % Thur II
RegioCovID(reach_upstream(reach_upstream(:,118)>0,118))=3; % Thur III
RegioCovID(reach_upstream(reach_upstream(:,132)>0,132))=4; % Thur IV
RegioCovID(reach_upstream(reach_upstream(:,141)>0,141))=5; % Thur V
RegioCovID(reach_upstream(reach_upstream(:,232)>0,232))=6; % Thur VI
RegioCovID(reach_upstream(reach_upstream(:,421)>0,421))=7; % Thur VII
RegioCovID(reach_upstream(reach_upstream(:,677)>0,677))=8; % Thur VIII
RegioCovID(reach_upstream(reach_upstream(:,322)>0,322))=9; % Glatt I
RegioCovID(reach_upstream(reach_upstream(:,297)>0,297))=9; % Glatt I
RegioCovID(reach_upstream(reach_upstream(:,476)>0,476))=10; % Glatt II
RegioCovID(reach_upstream(reach_upstream(:,499)>0,499))=11; % Wissbach
RegioCovID(reach_upstream(reach_upstream(:,163)>0,163))=12; % Necker I
RegioCovID(reach_upstream(reach_upstream(:,264)>0,264))=13; % Necker II
RegioCovID(reach_upstream(reach_upstream(:,396)>0,396))=14; % Necker III
RegioCovID(reach_upstream(reach_upstream(:,92)>0,92))=15; % Gonzenbach          
RegioCovID(reach_upstream(reach_upstream(:,113)>0,113))=16; % Dietfurterbach
RegioCovID(reach_upstream(reach_upstream(:,479)>0,479))=17; % Luteren
    
% build covariate matrix
RegioCovMat=zeros(N_reach,18);
for i=1:17
    RegioCovMat(RegioCovID==i,i)=1;
end
RegioCovMat(:,18)=AreaUpstream; % add area upstream
% Z normalize
RegioCovMat=(RegioCovMat-mean(RegioCovMat))./std(RegioCovMat);
