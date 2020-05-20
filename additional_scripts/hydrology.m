
% load network data, discharge duration curves, stage-discharge relationships
%load('ThurData.mat')
[DQ2181,~,~]=xlsread('StageDischarge.xlsx','2181'); % not used, it's downstream of the catchment outlet
[DQ2303,~,~]=xlsread('StageDischarge.xlsx','2303');
[DQ2374,~,~]=xlsread('StageDischarge.xlsx','2374');
[DQ2414,~,~]=xlsread('StageDischarge.xlsx','2414');
[DQ2305,~,~]=xlsread('StageDischarge.xlsx','2305');

[Qjun2016,~,~]=xlsread('ThurQ_jun2016.xlsx'); % from https://www.hydrodaten.admin.ch/en/2303.html (2374, 2414, 2305)
Qjun2016=Qjun2016(15,:);

widths=[40 35 20 2.5 8]; % from aerial images
Area=[1085 493 88.1 3.19 16.7]; % from https://www.hydrodaten.admin.ch/en/2181.html (2303, 2374, 2414, 2305)

% site coordinates (from FOEN data)
coord=[733560 263180; % 2181
       723675 252720; % 2303
       727110 247290; % 2374
       718840 248440; % 2414
       737270 251290]; % 2305

% estimate depth
depth=zeros(5,1);
depth(2)=interp1(DQ2303(:,2),DQ2303(:,1)-DQ2303(1,1),Qjun2016(1));
depth(3)=interp1(DQ2374(:,2),DQ2374(:,1)-DQ2374(1,1),Qjun2016(2));
depth(4)=interp1(DQ2414(:,2),DQ2414(:,1)-DQ2414(1,1),Qjun2016(3));
depth(5)=interp1(DQ2305(:,2),DQ2305(:,1)-DQ2305(1,1),Qjun2016(4));

% fit power-law relationships for Q,D,W, among sites 2023, 2374, 2414, 2305
% (exclude 2414 for D)
lmQ=fitlm(log(Area(2:5)),log(Qjun2016));
lmW=fitlm(log(Area(2:5)),log(widths(2:5)));
lmD=fitlm(log(Area([2 3 5])),log(depth([2 3 5])));

% derive power-law relationship for V
coeffV=exp(lmQ.Coefficients{1,1}-lmD.Coefficients{1,1}-lmW.Coefficients{1,1});
expV=lmQ.Coefficients{2,1}-lmD.Coefficients{2,1}-lmW.Coefficients{2,1};

%% figure power laws
color=get(gca,'ColorOrder'); close all;
figure('Units','Centimeters','Position',[0 0 30 8])
subplot(1,4,1)
loglog(Area(2:5),widths(2:5),'ok','markerfacecolor','k'); hold on;
loglog([0.1 1000],exp(lmW.Coefficients{1,1})*[0.1 1000].^(lmW.Coefficients{2,1}))
text(0.15,70,['R^2 = ',num2str(round(10000*lmW.Rsquared.Ordinary)/10000)],'color',color(1,:))
text(0.15,40,['w = ',num2str(round(1e3*exp(lmW.Coefficients{1,1}))/1e3),...
    'A^{',num2str(round(1e3*lmW.Coefficients{2,1})/1e3),'}'],'color',color(1,:))
plot([max(AreaUpstream) max(AreaUpstream)],[0.01 100],'--k')
plot([min(AreaUpstream) min(AreaUpstream)],[0.01 100],'--k')
xlabel('Contributing Area [km^2]'); ylabel('River width [m]')
set(gca,'tickdir','out','xlim',[0.1 1000],'xtick',[0.1 1 10 100 1000])
subplot(1,4,2)
loglog(Area(2:5),Qjun2016,'ok','markerfacecolor','k'); hold on;
loglog([0.1 1000],exp(lmQ.Coefficients{1,1})*[0.1 1000].^(lmQ.Coefficients{2,1}))
text(0.15,700,['R^2 = ',num2str(round(10000*lmQ.Rsquared.Ordinary)/10000)],'color',color(1,:))
text(0.15,300,['Q = ',num2str(round(1e3*exp(lmQ.Coefficients{1,1}))/1e3),...
    'A^{',num2str(round(1e3*lmQ.Coefficients{2,1})/1e3),'}'],'color',color(1,:))
plot([max(AreaUpstream) max(AreaUpstream)],[0.001 1000],'--k')
plot([min(AreaUpstream) min(AreaUpstream)],[0.001 1000],'--k')
xlabel('Contributing Area [km^2]'); ylabel('Water discharge [m^3/s]')
set(gca,'tickdir','out','xlim',[0.1 1000],'xtick',[0.1 1 10 100 1000])
subplot(1,4,3)
loglog(Area([2 3 5]),depth([2 3 5]),'ok','markerfacecolor','k'); hold on;
loglog([0.1 1000],exp(lmD.Coefficients{1,1})*[0.1 1000].^(lmD.Coefficients{2,1}))
text(0.15,8,['R^2 = ',num2str(round(10000*lmD.Rsquared.Ordinary)/10000)],'color',color(1,:))
text(0.15,5,['D = ',num2str(round(1e3*exp(lmD.Coefficients{1,1}))/1e3),...
    'A^{',num2str(round(1e3*lmD.Coefficients{2,1})/1e3),'}'],'color',color(1,:))
plot([max(AreaUpstream) max(AreaUpstream)],[0.01 100],'--k')
plot([min(AreaUpstream) min(AreaUpstream)],[0.01 100],'--k')
xlabel('Contributing Area [km^2]'); ylabel('Water depth [m]')
set(gca,'tickdir','out','xlim',[0.1 1000],'ylim',[0.01 10],'xtick',[0.1 1 10 100 1000])
subplot(1,4,4)
loglog([0.1 1000],coeffV*[0.1 1000].^expV); hold on
plot([max(AreaUpstream) max(AreaUpstream)],[0.01 100],'--k')
plot([min(AreaUpstream) min(AreaUpstream)],[0.01 100],'--k')
text(0.15,8,['v = ',num2str(round(1e3*coeffV)/1e3),...
    'A^{',num2str(round(1e3*expV)/1e3),'}'],'color',color(1,:))
xlabel('Contributing Area [km^2]'); ylabel('Water velocity [m/s]')
set(gca,'tickdir','out','xlim',[0.1 1000],'xtick',[0.1 1 10 100 1000],'ylim',[0.1 10])

% export csv
T=table(Area(2:5)',Qjun2016',widths(2:5)',[depth(2:3); NaN; depth(5)],...
    'VariableNames',{'A','Q','w','D'});
writetable(T,'source_data/FigS2.csv')

%%
% apply power-law relationships across the catchment
Qjun=exp(lmQ.Coefficients{1,1})*AreaUpstream.^(lmQ.Coefficients{2,1});
ReachWidth=exp(lmW.Coefficients{1,1})*AreaUpstream.^(lmW.Coefficients{2,1});
ReachDepth=exp(lmD.Coefficients{1,1})*AreaUpstream.^(lmD.Coefficients{2,1});
Velocity=coeffV*AreaUpstream.^expV;

save('data/ThurHydrology.mat','Qjun','ReachWidth','ReachDepth','Velocity')

