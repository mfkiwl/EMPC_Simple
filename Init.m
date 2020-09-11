clear all; clc; close all
%% adding some toolboxes
addpath C:\Werk_Data\Software\export_fig;
addpath C:\Werk_Data\Software\hline_vline;
addpath('C:\Program Files\Mosek\9.2\toolbox\R2015a\')
setenv('PATH', [getenv('PATH') ';C:\Program Files\Mosek\9.2\tools\platform\win64x86\bin']);
addpath Functions
global Parameters
%% loading linearization data
load 'NREL5MW_CPdata_new3'
Plot_Curves=0;
load Parameters

%% parameters of the turbine
Ts=0.5; Parameters.Ts=Ts;
Parameters.R=63;
Parameters.rho=1.225;
Parameters.G=97;
Parameters.J= 4.2e7;% 3.5444E+07
Parameters.Prated=5e6;
Parameters.Effic=0.944;
Parameters.Omegagmax=12.1/60*2*pi*Parameters.G*1.1; % allow 10% overspeeding
Parameters.Omegagrated=12.1/60*2*pi*Parameters.G; % allow 10% overspeeding
Parameters.Omegagmin=4.1/60*2*pi*Parameters.G;
Parameters.Tmax=Parameters.Prated/Parameters.Omegagmax/Parameters.Effic;
Parameters.ScaleP=Parameters.Prated;
Parameters.ScaleK=1/2*Parameters.J*Parameters.Omegagmax.^2;
s=tf('s'); omega_f=1*2*pi; filter=omega_f^2/(s^2+2*0.7*omega_f*s+omega_f^2); filter_d=c2d(filter,Ts);
%% creating some interpolation tables
Rotor_Lambda=Rotor_Lamda(1:150);
[Xq,Yq] = meshgrid(Rotor_Lambda, Rotor_Pitch);
Cp = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cP,Xq,Yq,'spline');Parameters.Rotor_cP=Rotor_cP; 
Cq = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cQ,Xq,Yq,'spline');Parameters.Rotor_cQ=Rotor_cQ; Parameters.Rotor_Lamda=Rotor_Lamda; Parameters.Rotor_Pitch=Rotor_Pitch;
Ct = interp2(Rotor_Lamda,Rotor_Pitch,Rotor_cT,Xq,Yq,'spline');
if Plot_Curves==1; Plot_Cp_Cq; end

%% compute optimal mode gain for Torque=k*omega^2
[II,JJ]=max(Cp);
[maxCp,III]=max(II); Pitch_opt=Yq(JJ(III),III); TSR_opt=Xq(JJ(III),III);
Kopt=1/2*Parameters.rho*pi*Parameters.R^5*maxCp/Parameters.G^3/TSR_opt^3/Parameters.Effic;
    
omega_0=0.05*2*pi; % init condition for rotor speed


%% form Constraints on PW
Parameters.Nconst=20; % number of lines

if Plot_Curves==1; 
    gg=0;
for V=5:1:25;
    gg=gg+1;
K=1/2*(Rotor_Lamda*V/Parameters.R*Parameters.G).^2*Parameters.J./Parameters.ScaleK;
Rotor_K=1/2*(Rotor_Lamda*V/Parameters.R*Parameters.G).^2*Parameters.J./Parameters.ScaleK;
Rotor_Kls=0:min(max(Rotor_K)*0.7,1/2*Parameters.J*Parameters.Omegagmax.^2*10)/2000:min(max(Rotor_K)*0.7,1/2*Parameters.J*Parameters.Omegagmax.^2*10);
Cpp = interp2(Rotor_K,Rotor_Pitch,Rotor_cP,Rotor_Kls,Pitch_opt,'spline');
P=1/2*Parameters.rho*pi*Parameters.R^2*Cpp*V^3/Parameters.ScaleP;
Pf=1/2*Parameters.rho*pi*Parameters.R^2*Cpp/Parameters.ScaleP;
figure(101)
plot(Rotor_Kls,Pf); hold on
for real=1:1:1000
  [ijj,ijjj]=max(Pf);
  QQ=[0.8:0.05:1.2];
FF=sort([(sort(abs(rand(Parameters.Nconst-length(QQ),1)))*max(Rotor_Kls));Rotor_Kls(ijjj)*QQ']);
Pff = interp1(Rotor_Kls,Pf,FF,'pchip');
Bc(:,real)=interp1(Rotor_Kls,gradient(Pf)./mean(diff(Rotor_Kls)),FF,'pchip');
Ac(:,real)=Pff-Bc(:,real).*FF;
%plot(FF,Pff,'+');
for ij=1:1:length(FF)
%plot(Rotor_Kls,Ac(ij)+Bc(ij)*Rotor_Kls,'k:')
Complete_set(:,ij,real)=Ac(ij,real)+Bc(ij,real)*Rotor_Kls;
end
%plot(Rotor_Kls,min(squeeze(Complete_set(:,:,real))')','k','Linewidth',2)
NORM(real)=vaf(min(squeeze(Complete_set(:,:,real))'),Pf);

end
[jj,jjj]=max(NORM); display(num2str(jj))
plot(Rotor_Kls,min(squeeze(Complete_set(:,:,jjj))')','k','Linewidth',2)
Parameters.BOUNDPW.V(gg)=V;
Parameters.BOUNDPW.Bc(:,gg)=Bc(:,jjj);
Parameters.BOUNDPW.Ac(:,gg)=Ac(:,jjj);
clear Complete_set;
plot(Rotor_Kls,sqrt(2*Rotor_Kls*Parameters.ScaleK/Parameters.J)*Parameters.Tmax/Parameters.ScaleP/V^3,'r','linewidth',3)
%plot(FF,Parameters.BOUNDPW.Bc(:,gg).*FF+Parameters.BOUNDPW.Ac(:,gg),'+')
end
save Parameters 'Parameters'
end
