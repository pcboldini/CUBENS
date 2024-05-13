close all
clear all
clc

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26);

% What you want 

Re_z1=550; % start location turbulent region
Re_z2=670; % end location turbulent region
pert_zLen=10; % start location of strip
pert_Remid=415; % end location of strip

zplus_min=10; % zplus turbulent region (DNS)
zpluspert_min=5; % zplus strip region (DNS)
zplus_max=20; % zplus laminar region (DNS)

nz=1024; % number of points
lz_start=0;
lz=250.3919; % length of domain
DNS_Re_99Sta=1557.55;
Restart=316.23;
ReTau=50.0175;

% Factors to adapt -> DNS
z_pert1fac=0.17; % start of strip (DNS)
z_pert2fac=0.17; % end of strip (DNS)
z_1fac=0.42; % start of turbulent (DNS)
z_2fac=1.2; % end of turbulent (DNS)

bumpzpert2=0.01; % bump end of strip (DNS)
bumpzpert1=0.01; % bump start of strip (DNS)
bumpz1=0.1; % bump start of turbulent (DNS)
bumpz2=0.04; % bump end of turbulent (DNS)


%% Calculation

deltapert1=z_pert1fac*bumpzpert1;
deltapert2=z_pert2fac*bumpzpert2;
delta1=z_1fac*bumpz1;
delta2=z_2fac*bumpz2;

for i=1:nz
    z_old(i)=(i-1)/(nz-1)*lz;
end

Deltaz_old=(z_old(2)-z_old(1))*ones(nz,1);

for i=1:nz
    index(i)=i;
    fact(i)   =  (i-1.0)/(nz-1);
    dz(i)=(zplus_min+0.5*(zplus_max-zplus_min)*( 2-tanh((fact(i)-z_1fac)/(delta1))+tanh((fact(i)-z_2fac)/(delta2)) )+0.5*(zplus_max-zpluspert_min)*(-tanh((fact(i)-z_pert1fac)/(deltapert1))+tanh((fact(i)-z_pert2fac)/(deltapert2)) ) );    
    %dz2(i)=zplus_min+0.5*(zplus_max-zplus_min)*( 1-tanh((fact(i)-z_1fac)/(delta1)) );    

    ddz(i)=0.5*(zplus_max-zplus_min)*( (sech((z_2fac-fact(i))/delta2)).^2/delta2 - (sech((z_1fac-fact(i))/delta1)).^2/delta1 )  + 0.5*(zplus_max-zpluspert_min)*( (sech((z_pert2fac-fact(i))/deltapert2)).^2/deltapert2 - (sech((z_pert1fac-fact(i))/deltapert1)).^2/deltapert1 ); % sech=1/cosh
end

z_int1=deltapert1*0.5*(zpluspert_min-zplus_max)*log(cosh((z_pert1fac)/deltapert1))...
        +deltapert2*0.5*(zplus_max-zpluspert_min)*log(cosh((z_pert2fac)/deltapert2))...
        +delta1*0.5*(zplus_min-zplus_max)*log(cosh((z_1fac)/delta1))...
        +delta2*0.5*(zplus_max-zplus_min)*log(cosh((z_2fac)/delta2));

z_intE=delta1*0.5*(zplus_min-zplus_max)*log(cosh((z_1fac-1)/delta1))...
      +delta2*0.5*(zplus_max-zplus_min)*log(cosh((z_2fac-1)/delta2))...
      +1*zplus_max;

z=zeros(nz,1);

for i=2:nz
    z(i)=deltapert1*0.5*(zpluspert_min-zplus_max)*log(cosh((z_pert1fac-fact(i))/deltapert1))...
        +deltapert2*0.5*(zplus_max-zpluspert_min)*log(cosh((z_pert2fac-fact(i))/deltapert2))...   
        +delta1*0.5*(zplus_min-zplus_max)*log(cosh((z_1fac-fact(i))/delta1))...
        +delta2*0.5*(zplus_max-zplus_min)*log(cosh((z_2fac-fact(i))/delta2))...
        +fact(i)*zplus_max-z_int1;
%     z11(i)= delta1*0.5*(zplus_min-zplus_max)*log(cosh((z_1fac-fact(i))/delta1))...
%          +fact(i)*zplus_max-z_int1;
end


scaling=z(end)/(lz);
z_new=z/scaling+lz_start; % scaling needed for the second derivative!
dz_new=dz/scaling; % scaling needed in the DNS
ddz_new=ddz/scaling; % scaling needed in the DNS

zplus=z_new*ReTau;

Deltaz=zeros(nz-1,1);
Delta_zplus=zeros(nz-1,1);

for i=1:nz-1
    Deltaz(i)=z_new(i+1)-z_new(i);
    Delta_zplus(i)=zplus(i+1)-zplus(i);
end

zstart=Restart^2/DNS_Re_99Sta;
z1=Re_z1^2/DNS_Re_99Sta-zstart;
z2=Re_z2^2/DNS_Re_99Sta-zstart;
[~,idx_z1] = min(abs(z_new-z1));
[~,idx_z2] = min(abs(z_new-z2));

pert_zmid=pert_Remid^2/DNS_Re_99Sta-zstart;
zpert1=pert_zmid-pert_zLen/2;
zpert2=pert_zmid+pert_zLen/2;

[~,idx_zpert1] = min(abs(z_new-zpert1));
[~,idx_zpert2] = min(abs(z_new-zpert2));

figure(1)
subplot(1,3,1)
plot(index,dz)
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ k $','interpreter','latex');
ylabel('$ d z^+ $','interpreter','latex');
line([idx_zpert1 idx_zpert1],[min(dz) max(dz)],'Color','red','LineStyle','--');
line([idx_zpert2 idx_zpert2],[min(dz) max(dz)],'Color','red','LineStyle','-.');
line([idx_z1 idx_z1],[min(dz) max(dz)],'Color','magenta','LineStyle','--');
line([idx_z2 idx_z2],[min(dz) max(dz)],'Color','magenta','LineStyle','-.');
%xlim([5 10])
%ylim([0 y_coord_max_Cf])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

subplot(1,3,2)
plot(index,z)
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ k $','interpreter','latex');
ylabel('$ z^+ $','interpreter','latex');
line([idx_zpert1 idx_zpert1],[min(z) max(z)],'Color','red','LineStyle','--');
line([idx_zpert2 idx_zpert2],[min(z) max(z)],'Color','red','LineStyle','-.');
line([idx_z1 idx_z1],[min(z) max(z)],'Color','magenta','LineStyle','--');
line([idx_z2 idx_z2],[min(z) max(z)],'Color','magenta','LineStyle','-.');

%xlim([5 10])
%ylim([0 y_coord_max_Cf])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

subplot(1,3,3)
plot(index,ddz)
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ k $','interpreter','latex');
ylabel('$ ddz^+ $','interpreter','latex');
%xlim([5 10])
%ylim([0 y_coord_max_Cf])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

figure(2)

subplot(1,3,1)
plot(index,z_new,index,z_old)
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ k $','interpreter','latex');
ylabel('$ z $','interpreter','latex');
line([idx_zpert1 idx_zpert1],[min(z_new) max(z_new)],'Color','red','LineStyle','--');
line([idx_zpert2 idx_zpert2],[min(z_new) max(z_new)],'Color','red','LineStyle','-.');
line([idx_z1 idx_z1],[min(z_new) max(z_new)],'Color','magenta','LineStyle','--');
line([idx_z2 idx_z2],[min(z_new) max(z_new)],'Color','magenta','LineStyle','-.');

%xlim([5 10])
%ylim([0 y_coord_max_Cf])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

subplot(1,3,2)
plot(index(1:end-1),Deltaz,index,Deltaz_old)
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ k $','interpreter','latex');
ylabel('$ \Delta z $','interpreter','latex');
line([idx_zpert1 idx_zpert1],[min(Deltaz) max(Deltaz)],'Color','red','LineStyle','--');
line([idx_zpert2 idx_zpert2],[min(Deltaz) max(Deltaz)],'Color','red','LineStyle','-.');
line([idx_z1 idx_z1],[min(Deltaz) max(Deltaz)],'Color','magenta','LineStyle','--');
line([idx_z2 idx_z2],[min(Deltaz) max(Deltaz)],'Color','magenta','LineStyle','-.');

%xlim([5 10])
%ylim([0 y_coord_max_Cf])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

subplot(1,3,3)
plot(index,ddz_new)
grid off
ax=gca;
ax.TickLabelInterpreter='latex';
xlabel('$ k $','interpreter','latex');
ylabel('$ \Delta z $','interpreter','latex');
%xlim([5 10])
%ylim([0 y_coord_max_Cf])
set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');


