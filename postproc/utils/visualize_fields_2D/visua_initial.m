close all
clear all
clc

%% Parameters 

global visua_2D visua_1cut cut_Re

visua_2D='no'; % yes/no
visua_1cut='yes';  % yes/no

USE_EOS='VdW'; % IG/VdW
cut_Re='no';  % yes/no
z_cut=0; % streamwise coordinate for the cut
Re_cut=1500; % Re_\delta for the cut

path_new="../../planes";

%%

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26)

if strcmp(USE_EOS,'IG')
initBL=load('../../initBL_param.mat');
elseif strcmp(USE_EOS,'VdW')
initBL=load('../../initBL_param.mat');
end

Re_deltaSta=initBL.BL_PARA.Re_deltaSta;
Re_deltaEnd=initBL.BL_PARA.Re_deltaEnd;
DNS_Re_99Sta=initBL.DNS_PARA.DNS_Re_99Sta;
z_start=initBL.DNS_PARA.DNS_x_Sta;
Pr_inf=initBL.BL_PARA.Pr_inf;
Ec_inf=initBL.BL_PARA.Ec_inf;

fileID = fopen(strcat(path_new,'/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_new,'/y.bin'));
y = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_new,'/z.bin'));
z = fread(fileID,'double');
fclose(fileID);

imax = length(x);
jmax = length(y);
kmax = length(z);

z_max=max(z); 
z_min=min(z);

fprintf('Streamwise coordinate from z=%.f to z=%.3f \n',z_min,z_max);
fprintf('or Re_delta from %.3f to %.3f \n\n',Re_deltaSta,Re_deltaEnd);

var = {'p','r','t','w','u','v','mu','ka','e'};

for i=1:numel(var)
    fname{i} = sprintf(strcat(path_new,'/ypl.1.',var{i},'.0000000.bin'));
    fileID(i)=fopen(fname{i});
    data{i} = fread(fileID(i),imax*jmax*kmax,'double');
    data{i} = reshape(data{i},imax,jmax,kmax);
    fclose(fileID(i));
end

    p_yslice=data{1}(:,1,:);
    Rho_yslice=data{2}(:,1,:);
    T_yslice=data{3}(:,1,:);
    W_yslice=data{4}(:,1,:);
    U_yslice=data{5}(:,1,:);
    Mu_yslice=data{7}(:,1,:)*DNS_Re_99Sta;
    Ka_yslice=data{8}(:,1,:)*DNS_Re_99Sta*Pr_inf*Ec_inf;
    E_yslice=data{9}(:,1,:);
    
if strcmp(visua_2D,'yes') 
    [Z,X] = meshgrid(z,x);   
elseif strcmp(visua_1cut,'yes')
    
    if strcmp(cut_Re,'yes')
        z_cut=Re_cut^2/DNS_Re_99Sta-z_start;
        [~,index_z_cut]=min(abs(z_cut-z));
    
        z_cut=z(index_z_cut);
        Re_cut=sqrt((z_cut+z_start)*DNS_Re_99Sta);
        
        fprintf('Extract variables at z=%.3f at Re_delta=%.3f \n',z_cut,Re_cut);
        
    elseif strcmp(cut_Re,'no')       
        [~,index_z_cut]=min(abs(z_cut-z));
    
        z_cut=z(index_z_cut);
        Re_cut=sqrt((z_cut+z_start)*DNS_Re_99Sta);
        fprintf('Extract variables at z=%.3f at Re_delta=%.3f \n',z_cut,Re_cut);
        
    end
   
    p_cut=p_yslice(:,1,index_z_cut);
    Rho_cut=Rho_yslice(:,1,index_z_cut);
    T_cut=T_yslice(:,1,index_z_cut);
    W_cut=W_yslice(:,1,index_z_cut);
    U_cut=U_yslice(:,1,index_z_cut);
    Mu_cut=Mu_yslice(:,1,index_z_cut);
    Ka_cut=Ka_yslice(:,1,index_z_cut);
    E_cut=E_yslice(:,1,index_z_cut);
    
    diff_bl=W_cut-0.99; [~,index_bl]=min(abs(diff_bl));
    delta99=x(index_bl);
    
end
    
    %% Postprocessing
if strcmp(visua_2D,'yes') 
    
    p_max=max(max(p_yslice));
    p_min=min(min(p_yslice));
    Rho_max=max(max(Rho_yslice));
    Rho_min=min(min(Rho_yslice));
    T_max=max(max(T_yslice));
    T_min=min(min(T_yslice));
    W_max=max(max(W_yslice));
    W_min=min(min(W_yslice));
    U_max=max(max(U_yslice));
    U_min=min(min(U_yslice));
    Mu_max=max(max(Mu_yslice));
    Mu_min=min(min(Mu_yslice));
    Ka_max=max(max(Ka_yslice));
    Ka_min=min(min(Ka_yslice));
    
elseif strcmp(visua_1cut,'yes')

end    
    %% Plot
    
    if strcmp(visua_2D,'yes') 
        
        yLim_max=5;
        nsteps=200;

        figure(1)
        set(0, 'DefaultFigureRenderer', 'painters'); 
            subplot(2,2,1)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(Z,X,reshape(Rho_yslice,imax,kmax),'LevelList',linspace(Rho_min-0.1,Rho_max+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel("$x/\delta_{99,in}$")
        ylabel("$y/\delta_{99,in}$")
        ylim([0 yLim_max]);
        colormap(jet)
        cb.Label.String = '$\rho^*/\rho^*_{\infty}$';
        cb.Label.Interpreter = 'latex';

        subplot(2,2,2)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(Z,X,reshape(T_yslice,imax,kmax),'LevelList',linspace(T_min-0.1,T_max+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel("$x/\delta_{99,in}$")
        ylabel("$y/\delta_{99,in}$")
        ylim([0 yLim_max]);
        colormap(jet)
        cb.Label.String = '$T^*/T^*_{\infty}$';
        cb.Label.Interpreter = 'latex';

        subplot(2,2,3)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(Z,X,reshape(W_yslice,imax,kmax),'LevelList',linspace(0,W_max+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel("$x/\delta_{99,in}$")
        ylabel("$y/\delta_{99,in}$")
        ylim([0 yLim_max]);
        colormap(jet)
        cb.Label.String = '$u^*/u^*_{\infty}$';
        cb.Label.Interpreter = 'latex';

        subplot(2,2,4)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(Z,X,reshape(U_yslice,imax,kmax)*1e3,'LevelList',linspace(0,U_max*1e3+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel("$x/\delta_{99,in}$")
        ylabel("$y/\delta_{99,in}$")
        ylim([0 yLim_max]);
        colormap(jet)
        cb.Label.String = '$10^{3} \times v^*/u^*_{\infty} $';
        cb.Label.Interpreter = 'latex';
    
    elseif strcmp(visua_1cut,'yes')
        
        Ymax=delta99*1.5;
        
        subplot(2,3,1)
        plot(W_cut,x,'-','LineWidth',3,'Color','blue') 
        hold on
        line([0 1.05],[delta99 delta99],'LineStyle',':','LineWidth',2,'Color','blue')
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ u^*/u^*_\infty $','interpreter','latex');
        ylabel("$y/\delta_{99,in}$");
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');
        axis([0 1.05 0 Ymax])

        subplot(2,3,2)
        plot(T_cut,x,'-','LineWidth',3,'Color','blue')
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ T^*/T^*_\infty $','interpreter','latex');
        axis([0 max(T_cut)+1 0 Ymax])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

        subplot(2,3,3)
        plot(Rho_cut,x,'-','LineWidth',3,'Color','blue')
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ \rho^*/\rho^*_\infty $','interpreter','latex');
        axis([0 1.05 0 Ymax])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');
        
        subplot(2,3,4)
        plot(U_cut*1e3,x,'-','LineWidth',3,'Color','blue')
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ 10^3 \times v^*/u^*_\infty $','interpreter','latex');
        ylabel("$y/\delta_{99,in}$");
        ylim([0 Ymax])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

        subplot(2,3,5)
        plot(Mu_cut,x,'-','LineWidth',3,'Color','blue')
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ \mu^*/\mu^*_\infty $','interpreter','latex');
        axis([0 max(Mu_cut)+1 0 Ymax])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

        subplot(2,3,6)
        plot(Ka_cut,x,'-','LineWidth',3,'Color','blue')
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ \kappa^*/\kappa^*_\infty $','interpreter','latex');
        axis([0 max(Ka_cut)+1 0 Ymax])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');
  
    end


