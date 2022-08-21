
clear
d_list=[100:10:510]/1000;
Npart=length(d_list);
wsl=zeros(1,Npart);%ws list
Rel=wsl;%Re list
res_Rel=wsl;%res_Re list

rho_p=1800;%particle density (kg/m3)
rho_f=1025;%fluid density (kg/m3)
mu_f=0.00109;% seawater viscosity in Pa s (from: https://www.engineeringtoolbox.com/sea-water-properties-d_840.html for temperature at 20 degrees, which is the mean Med temperature (https://www.watertemperature.org/Mediterranean-Sea--Geo.html) )
Re_start=100;%initial guess Reynold number
g=9.81;
exp1=0.25;
exp2=0.08;
exp3=-5.05;
tol=1.e-10;

for k=1:Npart
d_p=d_list(k)*0.001;
Psi=1; % I consider clay particles as spheres 
    
Re_old=Re_start;
Cd_old=(24/Re_old)*(((1-Psi)/Re_old)+1)^exp1+(24/Re_old)*0.1806*(Re_old^0.6459)*Psi^(-(Re_old^exp2))+0.4251/(1+(6880.95/(Re_old*(Psi^exp3))));
wt_old=sqrt((4*g*d_p*(rho_p-rho_f))/(3*Cd_old*rho_f));
for i=1:1000
    Cd_new=(24/Re_old)*(((1-Psi)/Re_old)+1)^exp1+(24/Re_old)*0.1806*(Re_old^0.6459)*Psi^(-(Re_old^exp2))+0.4251/(1+(6880.95/(Re_old*(Psi^exp3))));
    wt_new=sqrt((4*g*d_p*(rho_p-rho_f))/(3*Cd_new*rho_f));
    Re_new=(rho_f*wt_new*d_p)/mu_f;
    res_Re=abs(real(Re_new-Re_old));
    %fprintf(formatSpec,i,res_Re);
    if(res_Re<=tol)
        break
    else
        Re_old=Re_new;
        Cd_old=Cd_new;
        wt_old=wt_new;
        continue               
    end
end
Re=Re_new;
Cd=Cd_new;
wt=wt_new*1000;
% formatSpec = '\n Re = %10.3f\n Cd = %10.3f\n wt (mm s^-1) = %10.4f \n';
% fprintf(formatSpec,Re,Cd,wt);

wsl(k)=real(wt); %in mm/s
Rel(k)=real(Re);
res_Rel(k)=res_Re;

end

wsl=wsl/1000*86400; %in m/d