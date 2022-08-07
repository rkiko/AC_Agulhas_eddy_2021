chdir('~/GIT/AC_Agulhas_eddy_2021/Scripts')
clear
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Trajectories_Cyclo_filt.mat')
clearvars -except Fields_Trajectories
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Coloc_BioArgo.mat')

traj_eddy1_lon=Cell_Traj_InterestCyclo{1,5};
traj_eddy1_lat=Cell_Traj_InterestCyclo{1,6};
radius_max_eddy1=Cell_Traj_InterestCyclo{1,8};
radius_out_eddy1=Cell_Traj_InterestCyclo{1,12};
datenum_eddy1=Cell_Traj_InterestCyclo{1,7};
vel_azimuth_eddy1=Cell_Traj_InterestCyclo{1,9};

traj_eddy3_lon=Cell_Traj_InterestCyclo{3,5};
traj_eddy3_lat=Cell_Traj_InterestCyclo{3,6};
radius_max_eddy3=Cell_Traj_InterestCyclo{3,8};
radius_out_eddy3=Cell_Traj_InterestCyclo{3,12};
datenum_eddy3=Cell_Traj_InterestCyclo{3,7};
vel_azimuth_eddy3=Cell_Traj_InterestCyclo{3,9};

clear xticklabels
xticklabels=cell(4,1);xticklabels(1)=cellstr('Datenum');xticklabels(2)=cellstr('Lon');xticklabels(3)=cellstr('Lat');
xticklabels(4)=cellstr('Rad_max');xticklabels(5)=cellstr('Rad_out');xticklabels(6)=cellstr('Vel_Azimuth');
M=[datenum_eddy1;traj_eddy1_lon;traj_eddy1_lat;radius_max_eddy1;radius_out_eddy1;vel_azimuth_eddy1]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv')

M=[datenum_eddy3;traj_eddy3_lon;traj_eddy3_lat;radius_max_eddy3;radius_out_eddy3;vel_azimuth_eddy3]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy3.csv')

clear xticklabels
xticklabels=cell(4,1);xticklabels(1)=cellstr('Distance_Centroid');xticklabels(2)=cellstr('Rad_max');
xticklabels(3)=cellstr('Rad_out');xticklabels(4)=cellstr('sel_insideEddy');xticklabels(6)=cellstr('Vel_Azimuth');
M=[Distance_Centroid,Rad_max,Rad_out,Distance_Centroid<Rad_max];
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_partialFile_an64m.csv')


clearvars -except Id_traj Id_eddy_coloc
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Contours_max_cyclo.mat')

date_vec=datevec(date_num);
lontmp=CEs_max{53,1,1};

%%%%%%%%%%%%%%%
%I manually identify the eddy 2 (the one merging with the main eddy around
%the beginning of August 2021) by looking for the eddy whos centroid on the
%6 august is close to a location identified manually of Remi snapshot
clear
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Trajectories_Cyclo_filt.mat')

dateref=datenum(2021,8,6);
lon_ref=6;
lat_ref=-34;
delta=1.5;

i=327;ct=0;
for i=1:length(Cyclonic_Trajectories)
datetmp=Cyclonic_Trajectories{i,7};
idx=find(datetmp==dateref,1);
if isempty(idx),continue,end
ct=ct+1;
lontmp=Cyclonic_Trajectories{i,5};lontmp=lontmp(idx);
lattmp=Cyclonic_Trajectories{i,6};lattmp=lattmp(idx);

if (lontmp>=(lon_ref-delta))&&(lontmp<(lon_ref+delta))&&(lattmp>=(lat_ref-delta))&&(lattmp<(lat_ref+delta))
    idx_eddy=i;
end

end   

clearvars -except idx_eddy Cyclonic_Trajectories Fields_Trajectories dateref lon_ref lat_ref
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Contours_max_cyclo.mat')
idx=find(date_num==dateref);
delta=1;
i=1;
for i=1:200
lontmp=CEs_max{i,idx,1};
lattmp=CEs_max{i,idx,2};
sel=(lontmp>=(lon_ref-delta))&(lontmp<(lon_ref+delta))&(lattmp>=(lat_ref-delta))&(lattmp<(lat_ref+delta));
if sum(sel)>0
    fprintf('i %d sum sel %d\n',i,sum(sel))
    hold on, plot(lontmp,lattmp,'b')
end

end

%idx_eddy=498;
traj_eddy2_lon=Cyclonic_Trajectories{idx_eddy,5};
traj_eddy2_lat=Cyclonic_Trajectories{idx_eddy,6};
radius_max_eddy2=Cyclonic_Trajectories{idx_eddy,8};
radius_out_eddy2=Cyclonic_Trajectories{idx_eddy,12};
datenum_eddy2=Cyclonic_Trajectories{idx_eddy,7};
vel_azimuth_eddy2=Cyclonic_Trajectories{idx_eddy,9};

clear xticklabels
xticklabels=cell(4,1);xticklabels(1)=cellstr('Datenum');xticklabels(2)=cellstr('Lon');xticklabels(3)=cellstr('Lat');
xticklabels(4)=cellstr('Rad_max');xticklabels(5)=cellstr('Rad_out');xticklabels(6)=cellstr('Vel_Azimuth');
M=[datenum_eddy2;traj_eddy2_lon;traj_eddy2_lat;radius_max_eddy2;radius_out_eddy2;vel_azimuth_eddy2]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy2.csv')

load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Coloc_BioArgo.mat','Cell_Traj_InterestCyclo')

for j=1:3
if j==2
id_detection=Cyclonic_Trajectories{idx_eddy,2};
else
id_detection=Cell_Traj_InterestCyclo{j,2};
end
i=1;lon_cont=nan(200,length(id_detection));lat_cont=nan(200,length(id_detection));maxl=nan(length(id_detection),1);
for i=1:length(id_detection)
[a,b]=ind2sub(size(id_cyclo),id_detection(i));
l=length(CEs_max{a,b,1});
lon_cont(1:l,i)= CEs_max{a,b,1};
lat_cont(1:l,i)= CEs_max{a,b,2};
maxl(i)=l;

end
maxl=max(maxl);
lon_cont=lon_cont(1:maxl,:);
lat_cont=lat_cont(1:maxl,:);

M=array2table(lon_cont);writetable(M,sprintf('~/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy%d_lon_an64m.csv',j))
M=array2table(lat_cont);writetable(M,sprintf('~/GIT/AC_Agulhas_eddy_2021/Data/an64/CEs_max_eddy%d_lat_an64m.csv',j))

end

figure(1),clf,plot(traj_eddy2_lon,traj_eddy2_lat,'b'),hold on,plot(traj_eddy1_lon,traj_eddy1_lat,'r'),plot(traj_eddy3_lon,traj_eddy3_lat,'m')
