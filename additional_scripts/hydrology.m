
% load network data, discharge duration curves, stage-discharge relationships
%load('ThurData.mat')
[Qnum,~,~]=xlsread('ThurQ.xlsx');
[DQ2181,~,~]=xlsread('StageDischarge.xlsx','2181');
[DQ2303,~,~]=xlsread('StageDischarge.xlsx','2303');
[DQ2374,~,~]=xlsread('StageDischarge.xlsx','2374');
[DQ2414,~,~]=xlsread('StageDischarge.xlsx','2414');
[DQ2305,~,~]=xlsread('StageDischarge.xlsx','2305');

widths=[40 35 20 2.5 8]; % from aerial images
Area=Qnum(26,2:6);

days=Qnum(1:24,1);
Q_2181=Qnum(1:24,2); Q_2303=Qnum(1:24,3); Q_2374=Qnum(1:24,4); Q_2414=Qnum(1:24,5); Q_2305=Qnum(1:24,6);
Q_all=[Q_2181 Q_2303 Q_2374 Q_2414 Q_2305];

% site coordinates (from FOEN data)
coord=[733560 263180; % 2181
       723675 252720; % 2303
       727110 247290; % 2374
       718840 248440; % 2414
       737270 251290]; % 2305
% find site location in the network 
site_loc=nan(5,1);  reach_site=nan(5,1);
for i=2:5
    dist=sqrt((coord(i,1)-X).^2+(coord(i,2)-Y).^2);
    site_loc(i)=find(dist==min(dist));
    reach_site(i)=reach(site_loc(i));
end

% estimate median depth
depth(1)=interp1(DQ2181(:,2),DQ2181(:,1)-DQ2181(1,1),Q_2181(13));
depth(2)=interp1(DQ2303(:,2),DQ2303(:,1)-DQ2303(1,1),Q_2303(13));
depth(3)=interp1(DQ2374(:,2),DQ2374(:,1)-DQ2374(1,1),Q_2374(13));
depth(4)=interp1(DQ2414(:,2),DQ2414(:,1)-DQ2414(1,1),Q_2414(13));
depth(5)=interp1(DQ2305(:,2),DQ2305(:,1)-DQ2305(1,1),Q_2305(13));

% fit power-law relationships for Q,D,W, among sites 2023, 2374, 2414, 2305
lmQ=fitlm(log(Area(2:5)),log(Q_all(13,2:5)));
lmD=fitlm(log(Area(2:5)),log(depth(2:5)));
lmW=fitlm(log(Area(2:5)),log(widths(2:5)));

% derive power-law relationship for V
coeffV=exp(lmQ.Coefficients{1,1}-lmD.Coefficients{1,1}-lmW.Coefficients{1,1});
expV=lmQ.Coefficients{2,1}-lmD.Coefficients{2,1}-lmW.Coefficients{2,1};

% apply power-law relationships across the catchment
Qmedian=exp(lmQ.Coefficients{1,1})*AreaUpstream.^(lmQ.Coefficients{2,1});
ReachWidth=exp(lmW.Coefficients{1,1})*AreaUpstream.^(lmW.Coefficients{2,1});
ReachMedianDepth=exp(lmD.Coefficients{1,1})*AreaUpstream.^(lmD.Coefficients{2,1});
VelocityMedian=coeffV*AreaUpstream.^expV;

save('data/ThurHydrology.mat','Qmedian','ReachWidth','ReachMedianDepth','VelocityMedian')

