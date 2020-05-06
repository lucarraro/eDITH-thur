function [f1]=DrawRiverMap(vector,maxval,minval,title_string,cmap,geometry,newfig,show_leg,colmapZissou)

% load geometrical features of the river
[X,Y,XX,YY,Xc,Yc,subcatch,N_reach,AD_pixel,nnodes,outlet,area_upstream,reach]=v2struct(geometry);

% define colormaps
if strcmp(cmap,'FRA')
    colmap=[[linspace(0,0.9,50)'; linspace(0.9,1,50)'] [linspace(0,0.9,50)';
        linspace(0.9,0,50)'] [linspace(1,0.9,50)'; linspace(0.9,0,50)']];
elseif strcmp(cmap,'ITA')
    colmap=[[linspace(1,0.9,50)'; linspace(0.9,0,50)'] [linspace(0,0.9,50)';
        linspace(0.9,0.7,50)'] [linspace(0,0.9,50)'; linspace(0.9,0,50)']];
elseif strcmp(cmap,'GER')
    colmap=[[linspace(1,1,50)'; linspace(1,0,50)'] [linspace(1,0,50)';
        linspace(0,0,50)'] [linspace(0,0,50)'; linspace(0,0,50)']];
elseif strcmp(cmap,'JET')
    colmap=jet(100);
elseif strcmp(cmap,'B&W')
    colmap=[linspace(0,1,100)' linspace(0,1,100)' linspace(0,1,100)'];
elseif strcmp(cmap,'W&R')
    colmap=[linspace(1,1,100)' linspace(1,0,100)' linspace(1,0,100)'];
elseif strcmp(cmap,'W&B')
    colmap=[linspace(1,0,100)' linspace(1,0,100)' linspace(1,0,100)'];
elseif strcmp(cmap,'W&DR')
    colmap=[[linspace(1,1,50)'; linspace(1,0.5,50)'] [(linspace(1,0,50)).^0.5'; 
        linspace(0,0,50)'] [(linspace(1,0,50)).^0.5'; linspace(0,0,50)']];
elseif strcmp(cmap,'ZIS')
    ind=round(linspace(1,length(colmapZissou),100));
    colmap=colmapZissou(ind,:);
end

% create new figure
if newfig==2
    f1=figure('units','centimeters','position',[0 5 10 14],'visible','off');
elseif newfig==1
    f1=figure('units','centimeters','position',[0 5 10 14],'visible','on');
end
axis equal; axis off;
hold on; 

% draw catchment contour
for i=1:length(Xc)-1
    line([Xc(i) Xc(i+1)],[Yc(i) Yc(i+1)],'color',[0.4 0.4 0.4],'linewidth',0.5);
end
line([Xc(end) Xc(1)],[Yc(end) Yc(1)],'color',[0.4 0.4 0.4],'linewidth',0.5);

% draw map by filling subcatchments (needs definition of subcatchment borders Xb, Yb)
% if strcmp(type,'subc') % draw fig with subcatchments
%     for i=1:N_reach
%         if isnan(vector(i))==0
%             index=floor((vector(i)-minval)/(maxval-minval)*100);
%             if index==0; index=1; elseif index>100; index=100; end
%             %plot(XX(find(subcatch==i)),YY(find(subcatch==i)),'.','color',colmap(index,:));
%             patch(Xb.(['r',num2str(i)]),Yb.(['r',num2str(i)]),colmap(index,:),'edgecolor','none')
%         end
%     end
%     colormap(colmap); caxis([minval maxval]); title(title_string);
%     if show_leg; colorbar; end
%     if show_net
%         for i=1:nnodes
%             if i~=outlet
%                 ind_s=reach(i);
%                 line([X(i) X(AD_pixel(i,:)==1)],[Y(i) Y(AD_pixel(i,:)==1)],...
%                     'color','b','linewidth',0.9*(area_upstream(ind_s)*1e-1)^0.4,'AlignVertexCenters','on')
%             end
%         end       
%     end
% elseif strcmp(type,'netw') % draw fig with stretches
    
    colormap(colmap); caxis([minval maxval]); title(title_string);
    if show_leg; colorbar; end
    for i=1:nnodes
        if i~=outlet
            ind_s=reach(i);
            if isnan(vector(ind_s))==0
                indcol=floor((vector(ind_s)-minval)/(maxval-minval)*100);
                if indcol<=0; indcol=1; end 
                if indcol>100; indcol=100; end                
                line([X(i) X(AD_pixel(i,:)==1)],[Y(i) Y(AD_pixel(i,:)==1)],...
                    'color',colmap(indcol,:),'linewidth',0.9*(area_upstream(ind_s)*1e-1)^0.4,'AlignVertexCenters','on')
            else
                line([X(i) X(AD_pixel(i,:)==1)],[Y(i) Y(AD_pixel(i,:)==1)],...
                    'color','b','linewidth',0.9*(area_upstream(ind_s)*1e-1)^0.4,'AlignVertexCenters','on')
            end
        end
end


