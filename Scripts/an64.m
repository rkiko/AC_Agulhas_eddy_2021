chdir('~/GIT/AC_Agulhas_eddy_2021/Scripts')
clear
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Trajectories_Cyclo_filt.mat')
clear
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Coloc_BioArgo.mat')

traj_eddy1_lon=Cell_Traj_InterestCyclo{1,5};
traj_eddy1_lat=Cell_Traj_InterestCyclo{1,6};
radius_max_eddy1=Cell_Traj_InterestCyclo{1,8};
radius_out_eddy1=Cell_Traj_InterestCyclo{1,12};
datenum_eddy1=Cell_Traj_InterestCyclo{1,7};

traj_eddy3_lon=Cell_Traj_InterestCyclo{3,5};
traj_eddy3_lat=Cell_Traj_InterestCyclo{3,6};
radius_max_eddy3=Cell_Traj_InterestCyclo{3,8};
radius_out_eddy3=Cell_Traj_InterestCyclo{3,12};
datenum_eddy3=Cell_Traj_InterestCyclo{3,7};

clear xticklabels
xticklabels=cell(4,1);xticklabels(1)=cellstr('Datenum');xticklabels(2)=cellstr('Lon');xticklabels(3)=cellstr('Lat');
xticklabels(4)=cellstr('Rad_max');xticklabels(5)=cellstr('Rad_out');
M=[datenum_eddy1;traj_eddy1_lon;traj_eddy1_lat;radius_max_eddy1;radius_out_eddy1]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv')

M=[datenum_eddy3;traj_eddy3_lon;traj_eddy3_lat;radius_max_eddy3;radius_out_eddy3]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy3.csv')

clear xticklabels
xticklabels=cell(4,1);xticklabels(1)=cellstr('Distance_Centroid');xticklabels(2)=cellstr('Rad_max');
xticklabels(3)=cellstr('Rad_out');xticklabels(4)=cellstr('sel_insideEddy');
M=[Distance_Centroid,Rad_max,Rad_out,Distance_Centroid<Rad_max];
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/Distance_and_Radius_partialFile_an64m.csv')


clearvars -except Id_traj Id_eddy_coloc
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Contours_max_cyclo.mat')

date_vec=datevec(date_num);
lontmp=CEs_max{53,1,1};
