%% Define colormap for DTM
row_darkgreen=1; row_green=50; row_yellow=100; row_lightbrown=250; row_brown=400; row_white1=950; row_white2=1000;
cmapDTM=zeros(1000,3);
cmapDTM(row_darkgreen,:)=[0 0.4 0]; cmapDTM(row_green,:)=[0.3 0.6 0]; 
for i=row_darkgreen+1:row_green-1
    slope=(i-row_darkgreen)/(row_green-row_darkgreen);
    for j=1:3; cmapDTM(i,j)=cmapDTM(row_darkgreen,j)+(cmapDTM(row_green,j)-cmapDTM(row_darkgreen,j))*slope; end   
end
cmapDTM(row_yellow,:)=[0.9 0.9 0]; 
for i=row_green+1:row_yellow-1
    slope=(i-row_green)/(row_yellow-row_green);
    for j=1:3; cmapDTM(i,j)=cmapDTM(row_green,j)+(cmapDTM(row_yellow,j)-cmapDTM(row_green,j))*slope; end
end
cmapDTM(row_lightbrown,:)=[0.7 0.55 0]; 
for i=row_yellow+1:row_lightbrown-1
    slope=(i-row_yellow)/(row_lightbrown-row_yellow);
    for j=1:3; cmapDTM(i,j)=cmapDTM(row_yellow,j)+(cmapDTM(row_lightbrown,j)-cmapDTM(row_yellow,j))*slope; end
end
cmapDTM(row_brown,:)=[0.3 0.1 0]; 
for i=row_lightbrown+1:row_brown-1
    slope=(i-row_lightbrown)/(row_brown-row_lightbrown);
    for j=1:3; cmapDTM(i,j)=cmapDTM(row_lightbrown,j)+(cmapDTM(row_brown,j)-cmapDTM(row_lightbrown,j))*slope; end
end
cmapDTM(row_white1,:)=[0.9 0.9 0.9];
for i=row_brown+1:row_white1-1
    slope=(i-row_brown)/(row_white1-row_brown);
    for j=1:3; cmapDTM(i,j)=cmapDTM(row_brown,j)+(cmapDTM(row_white1,j)-cmapDTM(row_brown,j))*slope; end
end
cmapDTM(row_white2,:)=[0.9 0.9 0.9];
for i=row_white1+1:row_white2-1
    slope=(i-row_white1)/(row_white2-row_white1);
    for j=1:3; cmapDTM(i,j)=cmapDTM(row_white1,j)+(cmapDTM(row_white2,j)-cmapDTM(row_white1,j))*slope; end
end