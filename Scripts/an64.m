chdir('~/GIT/AC_Agulhas_eddy_2021/Scripts')
clear
load('~/GIT/AC_Agulhas_eddy_2021/Data/an64/Coloc_BioArgo.mat')

traj_eddy1_lon=Cell_Traj_InterestCyclo{1,5};
traj_eddy1_lat=Cell_Traj_InterestCyclo{1,6};

traj_eddy3_lon=Cell_Traj_InterestCyclo{3,5};
traj_eddy3_lat=Cell_Traj_InterestCyclo{3,6};

clear xticklabels
xticklabels=cell(2,1);xticklabels(1)=cellstr('Lon');xticklabels(2)=cellstr('Lat');

M=[traj_eddy1_lon;traj_eddy1_lat]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy1.csv')

M=[traj_eddy3_lon;traj_eddy3_lat]';
M=array2table(M,'VariableNames',xticklabels);writetable(M,'~/GIT/AC_Agulhas_eddy_2021/Data/an64/traj_eddy3.csv')


