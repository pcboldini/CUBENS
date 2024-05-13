clear all
clc
close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26)

%% Parameters 

prop = 'r'; % r / T / u / v / w / mu / ka / e 
cut='Re'; % Re/z/x
USE_EOS='IG'; % IG/VdW
path_new_planes="../../planes";
path_new_avg="../../../restart";

z_cut=20; % streamwise coordinate for the cut
Re_cut=1500; % Re_\delta for the cut
x_cut=1; % wall-normal coordinate for the cut

istart = 0000000;
iend = 0000100;
deltastep = 10;
compare_steps=[0];

%% Calculation

if strcmp(USE_EOS,'IG')
initBL=load('../../initBL_param.mat');
elseif strcmp(USE_EOS,'VdW')
initBL=load('../../initBL_param.mat');
end

Re_deltaSta=initBL.BL_PARA.Re_deltaSta;
Re_deltaEnd=initBL.BL_PARA.Re_deltaEnd;
DNS_Re_99Sta=initBL.DNS_PARA.DNS_Re_99Sta;
z_start=initBL.DNS_PARA.DNS_x_Sta;

fileID = fopen(strcat(path_new_planes,'/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_new_planes,'/z.bin'));
z = fread(fileID,'double');
fclose(fileID);
z = z - z(1);

if cut == 'Re'
     z_cut=Re_cut^2/DNS_Re_99Sta-z_start;
     [~,index_z_cut]=min(abs(z_cut-z));
    
     z_cut=z(index_z_cut);
     Re_cut=sqrt((z_cut+z_start)*DNS_Re_99Sta);
        
     fprintf('Extract variables at z=%.3f at Re_delta=%.3f \n',z_cut,Re_cut);        
elseif cut == 'z'       
     [~,index_z_cut]=min(abs(z_cut-z));
    
      z_cut=z(index_z_cut);
      Re_cut=sqrt((z_cut+z_start)*DNS_Re_99Sta);
      fprintf('Extract variables at z=%.3f at Re_delta=%.3f \n',z_cut,Re_cut);
      cut='Re';  
elseif cut=='x'
      [~,index_x_cut]=min(abs(x_cut-x));

      x_cut=x(index_x_cut);
      fprintf('Extract variables at x=%.3f\n',x_cut);

end

if cut == 'Re'
    BL_scaling=sqrt((z(index_z_cut)+initBL.DNS_PARA.DNS_x_Sta)./initBL.DNS_PARA.DNS_x_Sta);
    diff_bl=BL_scaling-x; [~,index_BL_scaling]=min(abs(diff_bl));
    x99=x(index_BL_scaling);
end

nx = length(x);
nz = length(z);

z_max=max(z); 
z_min=min(z);
x_max=max(x); 
x_min=min(x);

fprintf('Streamwise coordinate from z=%.f to z=%.3f \n',z_min,z_max);
fprintf('Wall-normal coordinate from x=%.f to x=%.3f \n',x_min,x_max);

A = zeros(nx,nz);
index=1;

fname = sprintf('../../../restart/Yave_%s.bin', prop);
fileID = fopen(fname);
mat_read = fread(fileID,nx*1*nz,'double');
fclose(fileID);

mat_read = reshape(mat_read,nx,1,nz);
    
A=mat_read;
       
if (cut=='x')
    A_plot=A(index_x_cut,:);
elseif cut=='Re'
    A_plot=A(:,index_z_cut);     
end

A_max=max(A_plot);
A_min=min(A_plot);
   
A_limit=[0.99*A_min 1.01*A_max];
    
switch cut
    case 'x'
        plot(z(1:nz),A_plot,'LineWidth',2)
        xlabel('$z/\delta_\mathrm{in}$')
        ylabel(prop)  
        title("$x/\delta_\mathrm{in}=$"+x_cut);
        set(gca,'YLim',A_limit,'Xlim',[min(z), max(z)])
        hold off
        pause(0.3)
    case 'Re'     
        plot(A_plot,x(1:nx),'LineWidth',2)
        ylabel('$x/\delta_\mathrm{in}$')
        xlabel(prop)
        title("$Re_\delta=$"+Re_cut);
        set(gca,'XLim',A_limit,'Ylim',[min(x), x99*2])
        hold off
        pause(0.3)
end
    