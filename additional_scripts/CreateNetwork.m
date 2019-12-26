
cmap_DTM; % load custom colormap for displaying DTM

%% Load data
% read DTM and header lines
fileID = fopen('DTM.asc');
DTM=textscan(fileID,'%f','HeaderLines',6); %f %s %s %f %f %f
frewind(fileID); % reset counter to 1
numLines = 6; % read 6 headerlines
text = cell(numLines,2);
 for i = 1:numLines
     text(i,:) = strsplit(fgetl(fileID)); 
     eval([text{i,1},'=str2double(text{i,2});']);
 end
fclose(fileID); clear text

% load DTMfilled
DTMfel=imread('DTMfel.tif'); DTMfel=double(DTMfel); DTMfel=flipud(DTMfel); 

% load FD (flow-direction matrix obtained with D8 algorithm)
FD=imread('DTMp.tif'); FD=double(FD); FD=flipud(FD); 

% load CA (contributing area matrix given outlet position and FD matrix)
CA=imread('DTMad8.tif'); CA=double(CA); CA=flipud(CA); 
CA(CA==-1)=NaN;

% change format and rotate DTM
DTM=cell2mat(DTM); DTM=reshape(DTM',ncols,nrows); DTM=fliplr(DTM); DTM=DTM';

% build coordinate matrices 
XX=repmat(xllcorner+[1:ncols]*cellsize,size(DTM,1),1);
YY=repmat((yllcorner+[1:nrows]*cellsize)',1,size(DTM,2));

% make contour
MASK=not(isnan(CA)); C=contourc(double(MASK),[1 1]);
Xc=XX(1,C(1,2:5890)'); Yc=YY(C(2,2:5890)',1)';  % adjust manually to remove double points in C

%% Stream definition by threshold
A_thr=800; % no. pixel identifying channelized area

RN=CA>A_thr;
index_net=find(RN==1); %index of the complete matrix that are network pixel 
[row col]=ind2sub(size(FD),index_net);

X=XX(index_net)-0.5*cellsize;
Y=YY(index_net)-0.5*cellsize;
Z=DTM(index_net); % ELEVATION FOR EVERY NETWORK NODE 
Zfel=DTMfel(index_net); 

Drainage_area=CA(index_net);
nnodes=length(X);

%% Build adjacency matrix
length_pixel=ones(nnodes,1)*cellsize;  %the length of a pixel can be equal to cellsize or cellsize*sqrt(2)
length_pixel(FD(index_net)==2 | FD(index_net)==4 | FD(index_net)==6 | FD(index_net)==8)=sqrt(2)*cellsize;

AD_pixel=sparse(nnodes,nnodes);
down_node=zeros(nnodes,1); check_p=[];
for p=1:nnodes 
    [mov_col,mov_row]=neigh(FD(index_net(p)));
    %d=find((abs(row-(row(p)+mov_row))<0.1 & abs(col-(col(p)+mov_col))<0.1));
    d=find(index_net==index_net(p)+mov_row+nrows*mov_col);
    if ~isempty(d)
        down_node(p)=d;
    end
    AD_pixel(p,d)=1;
    if Zfel(d)>Zfel(p)
        check_p=[check_p; p];
    end
end
degreeout=full(sum(AD_pixel,2));  
degreein=full(sum(AD_pixel));  
confluence=zeros(nnodes,1);
confluence(degreein>=2)=1;  %confluence==1 iff the node is a confluence, ==0 otherwise
source=zeros(nnodes,1);
source(degreein==0)=1;  %source==1 iff the node is a source, ==0 otherwise
outlet=find(degreeout==0);
length_pixel(outlet)=0; % why??

up_nodes=zeros(nnodes,max(degreein));
for p=1:nnodes  
    up_nodes(p,1:degreein(p))=find(AD_pixel(:,p));
end

% compute pixel slope
slope=zeros(nnodes,1);
slope(down_node>0)=((Zfel(down_node>0)-Zfel(down_node(down_node>0)))./length_pixel(down_node>0));

%% COMPUTE TOTAL UPSTREAM LENGTH, #CONFLUENCES UPSTREAM, STRAHLER ORDER, AVERAGE REACH SLOPE
%total upstrem length and distance to the outlet
total_upstream_length=zeros(nnodes,1);
d2out=zeros(nnodes,1);
n_confluences_upstream=zeros(nnodes,1);  %nuber of confluences upstream of each node

for node=1:nnodes
    %node
    pos=node;
    for iter=1:nnodes 
        if pos==outlet
            break
        end
        d2out(node)=d2out(node)+length_pixel(pos);
        pos=down_node(pos);
        total_upstream_length(pos)=total_upstream_length(pos)+length_pixel(node);
        n_confluences_upstream(pos)=n_confluences_upstream(pos)+confluence(node);
    end
end

%compute strahler order
[sort_d2out index_d2out]=sort(d2out,'descend');
stream_order=zeros(nnodes,1);
for cont=1:nnodes
    node=index_d2out(cont);
    if degreein(node)==0
        stream_order(node)=1;
    elseif degreein(node)==1
        stream_order(node)=stream_order(up_nodes(node,1));
    elseif degreein(node)==3
        stream_order(node)=2; %exception because there is one such node in this network
    else  %degreein(node)==2
        if stream_order(up_nodes(node,1))==stream_order(up_nodes(node,2))
            stream_order(node)=stream_order(up_nodes(node,1))+1;
        else
            stream_order(node)=max(stream_order(up_nodes(node,1:degreein(node))));
        end
    end
end

% define reaches
conf_out=confluence+double(degreeout==0); sour_conf=source+confluence;
end_node=zeros(nnodes,1); reach=zeros(nnodes,1);
k=1;
for i=1:nnodes
    if sour_conf(i)==1
        reach(i)=k; j=down_node(i); %*(source(i)==1);
        while conf_out(j)==0
           reach(j)=k; j=down_node(j); 
        end
        end_node(i)=j; k=k+1;
    end
end
reach(outlet)=reach(find(down_node==outlet));
N_reach=k-1;

%% identify subcatchments 
subcatch=nan(size(CA));
subcatch(index_net)=reach; % attribute reach ID to pixels belonging to the river network
for i=1:N_reach
    reach_pixels.(['r',num2str(i)])=index_net(find(reach==i));
end
ind_headpixels=find(CA==1);
for i=1:length(ind_headpixels)
    p=ind_headpixels(i);
    k=NaN; p_new=p; sub_p=[];
    while isnan(k)
        k=subcatch(p_new); % check if pixel has already been attributed a reach ID
        if isnan(k) 
            sub_p=[sub_p; p_new]; % if not, add to the current path
        [mov_col,mov_row]=neigh(FD(p_new)); % build path
        p_new=p_new+mov_row+nrows*mov_col;
        end
    end
    subcatch(sub_p)=k; % attribute reach ID to current path
    reach_pixels.(['r',num2str(k)])=[reach_pixels.(['r',num2str(k)]); sub_p];
end

%% reach length
length_reach=zeros(N_reach,1);  
for j=1:N_reach
    for i=1:nnodes
        if reach(i)==j
            length_reach(j)=length_reach(j)+length_pixel(i);
        end
    end
end

% find down_reach
down_reach=zeros(N_reach,1);
for j=1:nnodes
    if down_node(j)==0
        down_reach(reach(j))=0;
    elseif reach(down_node(j))~=reach(j)
        down_reach(reach(j))=reach(down_node(j));
    end
end
outlet_reach=find(down_reach==0);

% matrix containing reach IDs upstream of each reach
reach_upstream=zeros(N_reach); connect_in=zeros(3,N_reach);
for j=1:N_reach
    tmp=find(down_reach==j); % find 1st order upstream reaches
    if isempty(tmp)==0   
       reach_upstream(1,j)=tmp(1); 
       connect_in(1,j)=tmp(1);
       if numel(tmp)>1
          reach_upstream(2,j)=tmp(2);
          connect_in(2,j)=tmp(2);
          if numel(tmp)>2
          reach_upstream(3,j)=tmp(3);
          connect_in(3,j)=tmp(3);
          end
       end
    for i=1:N_reach % scroll down the column
        ind=find(reach_upstream(:,j)==0,1,'first'); % find first free slot to add further reaches
        tmp=find(down_reach==reach_upstream(i,j)); % find higher order upstream reaches
        if isempty(tmp)==0
            reach_upstream(ind,j)=tmp(1);
            if numel(tmp)>1; reach_upstream(ind+1,j)=tmp(2); 
                if numel(tmp)>2; reach_upstream(ind+2,j)=tmp(3); 
                end; end
        end
    end
    end
end
reach_upstream(reach_upstream==outlet_reach)=0; % correct because down_reach(outlet)=0
for j=1:N_reach
    reach_upstream(find(reach_upstream(:,j)==0,1,'first'),j)=j; % add the reach itself
end
n_reach_upstream=sum(reach_upstream>0);

%% reach_slope
reach_slope=zeros(N_reach,1);
for k=1:N_reach
    tmp=0;
    for i=1:nnodes
    if reach(i)==k
       tmp=tmp+slope(i)*length_pixel(i);
    end
    end
    reach_slope(k)=tmp/length_reach(k);
end
reach_slope(reach_slope==0)=5e-4; % attribute non-null slope to flat reaches

%% AD reach
AD_reach=zeros(N_reach);
for i=1:N_reach
    for k=1:N_reach
        if down_reach(i)==k
            AD_reach(i,k)=1;
        end
    end
end

%% length downstream
n_reach_downstream=zeros(N_reach); length_downstream=zeros(N_reach);
list_reach_downstream=cell(N_reach);
% 1st order reaches
for i=1:N_reach
    length_downstream(i,i)=length_reach(i);
    list_reach_downstream{i,i}=i;
    for j=1:N_reach
        if down_reach(i)==j
            n_reach_downstream(i,j)=1;
            length_downstream(i,j)=length_reach(i)+length_reach(j);
            list_reach_downstream{i,j}=[i; j];
        end     
end; end
% all other reaches
k=1;
while k<N_reach
for i=1:N_reach; for j=1:N_reach
    if n_reach_downstream(i,j)==k 
        tmp=find(n_reach_downstream(:,i)>0);
        if ~isempty(tmp)
            n_reach_downstream(tmp,j)=k+1;
            length_downstream(tmp,j)=length_downstream(i,j)+length_reach(tmp);
            for kkk=1:numel(tmp)
            list_reach_downstream{tmp(kkk),j}=[list_reach_downstream{i,j}; tmp(kkk)];
            end
        end
    end
end; end
k=k+1;
end

%% subcatchment statistics
LocalElevation=zeros(N_reach,1); AreaLocal=zeros(N_reach,1); AreaUpstream=zeros(N_reach,1); 
for i=1:N_reach
    LocalElevation(i)=mean(DTMfel(reach_pixels.(['r',num2str(i)]))); % m a.s.l.
    AreaLocal(i)=numel(reach_pixels.(['r',num2str(i)]))*cellsize^2/1e6; %km2
end
for i=1:N_reach
    AreaUpstream(i)=sum(AreaLocal(reach_upstream(reach_upstream(:,i)>0,i))); % km2
end

%% Stream order
[~,ind_downstream]=sort(n_reach_upstream);
stream_order_reach=zeros(N_reach,1);
for i=1:N_reach
    j=ind_downstream(i);
    tmp=find(down_reach==j);
    if isempty(tmp)
        stream_order_reach(j)=1;
    elseif length(find(stream_order_reach(tmp)==max(stream_order_reach(tmp))))>1
         stream_order_reach(j)=max(stream_order_reach(tmp))+1;
    else stream_order_reach(j)=max(stream_order_reach(tmp));
    end
end

%% Compare shapefile streams and sites
Sites=xlsread('SitesFunWorks_coordinates.xlsx');
siteID=Sites(:,1); siteX=Sites(:,2); siteY=Sites(:,3);

for i=1:numel(siteID)
    tentative_p(i)=find((sqrt((siteX(i)-X).^2+(siteY(i)-Y).^2))==min(sqrt((siteX(i)-X).^2+(siteY(i)-Y).^2)));
end

% manually correct site location (to prevent mismatches between site
% coordinates and position of the extracted river network)
% 22
j=tentative_p(7);k=1;
while reach(j)==reach(tentative_p(7))
    j=down_node(j);
    k=k+1;
end
tentative_p(7)=j;
%28
tentative_p(14)=7964; % checked on the map
% 30
j=tentative_p(5);k=1;
while reach(j)==reach(tentative_p(5))
    j=down_node(j);
    k=k+1;
end
tentative_p(5)=j;
%32
tentative_p(13)=14317; % checked on the map
%42
tentative_p(43)=15113; % checked on the map
%55
tentative_p(59)=22325; % checked on the map
%57
tentative_p(49)=15106; % checked on the map

SitesLocation=tentative_p; clear tentative_p
SitesReach=reach(SitesLocation);
      
%% Average upstream slope
slope_upstream=zeros(N_reach,1);
for i=1:N_reach
    sset=reach_upstream(reach_upstream(:,i)>0,i);
    slope_upstream(i)=sum(reach_slope(sset).*length_reach(sset))/sum(length_reach(sset));
end

%% subcatchment adjacency
SubAD=zeros(N_reach);
for i=1:N_reach
    tmp=reach_pixels.(['r',num2str(i)]); % pick all pixels from a subcatchment
    for k=1:length(tmp)
        ind=tmp(k);
        sset=[ind+nrows; ind-1; ind-nrows; ind+1]; % find neighbouring pixels
        neighbor_sub=subcatch(sset); % find to which subcatchment these pixels belong
        neighbor_sub(isnan(neighbor_sub))=[]; % remove neighbours outside catchment
        border=find(neighbor_sub~=i); % find neighbouring pixel at a border between subcatchments
        if ~isempty(border) 
            SubAD(i,neighbor_sub(border))=1; % add connection
        end
    end
end

%% Call Drawing function
geometry=v2struct(X,Y,XX,YY,Xc,Yc,subcatch,N_reach,AD_pixel,nnodes,outlet,AreaUpstream,reach);
CalibSites=1:length(SitesReach);

% Draw catchment map
DTMmask=DTM.*MASK; DTMmask(DTMmask==0)=NaN;
figure('units','centimeters','position',[0 0 30 30])
pcolor(XX,YY,DTMmask); shading flat; 
caxis([400 2400]); colormap(cmapDTM); 
axis equal; axis off; 
colorbar
hold on; 
for i=1:nnodes
    if i~=outlet
        line([X(i) X(AD_pixel(i,:)==1)],[Y(i) Y(AD_pixel(i,:)==1)],...
            'color','b','linewidth',0.9*(AreaUpstream(reach(i))*1e-1)^0.4,'AlignVertexCenters','on')
    end
end
plot(Xc,Yc,'k')
plot(X(SitesLocation(CalibSites)),Y(SitesLocation(CalibSites)),'or','markerfacecolor','r')

DischargeCoord=[733560 263180; % 2181
       723675 252720; % 2303
       727110 247290; % 2374
       718840 248440; % 2414
       737270 251290]; % 2305
   
plot(DischargeCoord(:,1),DischargeCoord(:,2),'^c','markerfacecolor','c')

save('data/ThurData.mat','AD_reach','AreaLocal','AreaUpstream','down_reach','length_downstream','length_reach','LocalElevation','list_reach_downstream',...
    'N_reach','reach_pixels','reach_slope','reach_upstream','siteID','SitesReach','siteX','siteY','slope_upstream','stream_order_reach',...
    'X','Y','XX','YY','Xc','Yc','subcatch','N_reach','AD_pixel','nnodes','outlet','reach')