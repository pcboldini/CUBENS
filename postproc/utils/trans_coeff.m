clear all
clc
% close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',24);

%% Import 

path_folder=''; 

initBL=load(strcat('../',path_folder,'/initBL_param.mat')); 
transCoeff=importdata(strcat('../',path_folder,'/wall_prop.txt')); 

fileID = fopen(strcat('../',path_folder,'/planes/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat('../',path_folder,'/planes/y.bin'));
y = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat('../',path_folder,'/planes/z.bin'));
z = fread(fileID,'double');
fclose(fileID);

imax = length(x);
jmax = length(y);
kmax = length(z);

%% Proc

Re_deltaSta=250;
Re_deltaEnd=850;
delta_eta_BL=initBL.BL_PARA.delta_eta_BL;
DNS_Re_99Sta=Re_deltaSta*delta_eta_BL;
z_start = Re_deltaSta^2/DNS_Re_99Sta;

Rho_inf=initBL.BL_PARA.rho_inf;
Mu_inf=initBL.BL_PARA.mu_inf;
Pr_inf=initBL.BL_PARA.Pr_inf;
Ec_inf=initBL.BL_PARA.Ec_inf;
gamma_inf=initBL.BL_PARA.gamma;
M_inf=initBL.BL_PARA.M_inf;
T_inf=initBL.BL_PARA.T_inf;
Cp=initBL.BL_PARA.cp;
Rhow=transCoeff.data(:,3);
Tw=transCoeff.data(:,2);
Muw=transCoeff.data(:,4);
Kaw=transCoeff.data(:,5);

Re_delta=sqrt((z+z_start).*DNS_Re_99Sta);
Re_xSta=Re_deltaSta.^2;
Re_xEnd=Re_deltaEnd.^2;
Re_x=Re_delta.^2;
rec=Pr_inf^(1/3);
T_aw=(1+rec*0.5*(gamma_inf-1)*M_inf^2);

deltaHw=Cp.*(T_aw-Tw);

FLOW=initBL.FLOW;
THERMO=initBL.THERMO;
y_DNS=initBL.BL_MESH.y_DNS;
Dy_BL=initBL.BL_MESH.Dy_BL;

[Re_BL,Cf_BL,St_BL] = extract_selfsim(Dy_BL,FLOW,THERMO,deltaHw,initBL.BL_PARA);

%% Extract data

fname = sprintf(strcat('../',path_folder,'/results/Yavg/Yave_T13.bin'));
fileID = fopen(fname);
data = fread(fileID,imax*kmax,'double');
fclose(fileID);
tau_xz = reshape(data,imax,kmax);

fname = sprintf(strcat('../',path_folder,'/results/Yavg/Yave_q1.bin'));
fileID = fopen(fname);
data = fread(fileID,imax*kmax,'double');
fclose(fileID);
q_x = reshape(data,imax,kmax);

Cf_DNS=tau_xz(1,:);
St_DNS=q_x(1,2:end)./deltaHw(2:end)';

%% Skin friction

Re_x_theo=linspace(1e3,5e6,200)';
C_w_theo=sqrt(Rhow(1).*Muw(1));
Cf_lam_theo = 0.664.*C_w_theo./sqrt(Re_x_theo);

a_turb=(rec*0.5*(gamma_inf-1)*M_inf^2./Tw(1)).^(0.5);
b_turb=(T_aw./Tw(1)-1);
A_turb=(2*a_turb.^2-b_turb)./(b_turb.^2+4*a_turb.^2).^(0.5);
B_turb=(b_turb)./(b_turb.^2+4*a_turb.^2).^(0.5);

func = @(x) (asin(A_turb)+asin(B_turb))./(x.*(T_aw-1)).^(0.5)-4.15*log10(Re_x_theo.*x./Muw(1))+1.7;
Cf_turb_theo =fsolve(func,Cf_lam_theo);

%% Heat transfer

St_lam_theo=0.5*Cf_lam_theo./Pr_inf.^(2/3);
St_turb_theo=0.5*Cf_turb_theo./Pr_inf.^(2/3);

%% Plot

x_coord_max=6.2;
y_coord_max=3.5;

figure(1)
subplot(1,2,1)
h1=plot(Re_BL.^2*1e-5,Cf_BL*1e3,'--','LineWidth',2,'Color','black');
hold on
h2=plot(Re_x_theo*1e-5,Cf_turb_theo*1e3,'-.','LineWidth',2,'Color','black');
h3=plot(Re_x(1:5:end)*1e-5,2*Cf_DNS(1:5:end)*1e3,'o','LineWidth',3,'Color','blue');
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ Re_x/10^5 $','interpreter','latex');
ylabel('$  C_f \times 10^{-3} $','interpreter','latex');
legend([h1 h2 h3],'Laminar','Turbulent','DNS','interpreter','latex');
xlim([1 x_coord_max])
ylim([0.4 y_coord_max])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

subplot(1,2,2)
h1=plot(Re_BL.^2*1e-5,St_BL*1e4,'--','LineWidth',2,'Color','black');
hold on
h2=plot(Re_x_theo*1e-5,St_turb_theo*1e4,'-.','LineWidth',2,'Color','black');
h3=plot(Re_x(2:end)*1e-5,St_DNS*1e4,'-','LineWidth',3,'Color','blue');
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ Re_x/10^5 $','interpreter','latex');
ylabel('$  St \times 10^{-3} $','interpreter','latex');
%legend([h1 h2 h3],'Laminar','Turbulent','DNS','interpreter','latex');
xlim([0 x_coord_max])
ylim([0 y_coord_max])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');





