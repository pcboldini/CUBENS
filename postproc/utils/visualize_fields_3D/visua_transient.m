close all
clear all
clc

%% Parameters 

global visua cut_Re cut_plane prop

prop='rho'; % r,u,v,w,e

visua='1D'; % 2D/1D
cut_plane='y';  % x/y/z
cut_coord=0; % x_cut/y_cut/z_cut

cut_Re='no';  % yes/no
Re_cut=1500; % Re_\delta for the cut

USE_EOS='IG'; % IG/VdW

path_new_planes="../../planes";
path_new_restart="../../../restart";

%%

fileStep  = 10;
fileStart = 0000000;
fileEnd   = 0000100-fileStep;

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

fileID = fopen(strcat(path_new_planes,'/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_new_planes,'/y.bin'));
y = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_new_planes,'/z.bin'));
z = fread(fileID,'double');
fclose(fileID);

imax = length(x);
jmax = length(y);
kmax = length(z);
ntot=imax*jmax*kmax;

z_max=max(z); 
z_min=min(z);

fprintf('Streamwise coordinate from z=%.f to z=%.3f \n',z_min,z_max);
fprintf('or Re_delta from %.3f to %.3f \n\n',Re_deltaSta,Re_deltaEnd);

if strcmp(visua,'2D') 
    if strcmp(cut_plane,'y')
        [Z,X] = meshgrid(z,x); 
        [~,index_cut]=min(abs(cut_coord-y));
        y_cut=y(index_cut);
        fprintf('Extract variables at y=%.3f \n',y_cut);
    elseif strcmp(cut_plane,'x')
        [Z,Y] = meshgrid(z,y); 
        [~,index_cut]=min(abs(cut_coord-x));
        x_cut=x(index_cut);
        fprintf('Extract variables at x=%.3f \n',x_cut);
    elseif strcmp(cut_plane,'z')
        [Y,X] = meshgrid(y,x); 
        [~,index_cut]=min(abs(cut_coord-z));
        z_cut=z(index_cut);
        fprintf('Extract variables at z=%.3f \n',z_cut);
    end

elseif strcmp(visua,'1D') 
    if strcmp(cut_plane,'z')
        if strcmp(cut_Re,'yes')
            z_cut=Re_cut^2/DNS_Re_99Sta-z_start;
            [~,index_z_cut]=min(abs(z_cut-z));
    
            z_cut=z(index_z_cut);
            Re_cut=sqrt((z_cut+z_start)*DNS_Re_99Sta);
        
            fprintf('Extract variables at z=%.3f at Re_delta=%.3f \n',z_cut,Re_cut);
        
        elseif strcmp(cut_Re,'no')       
            [~,index_z_cut]=min(abs(cut_coord-z));
            z_cut=z(index_z_cut);
            Re_cut=sqrt((z_cut+z_start)*DNS_Re_99Sta);
            fprintf('Extract variables at z=%.3f at Re_delta=%.3f \n',z_cut,Re_cut);
        
        end

    elseif strcmp(cut_plane,'x')

      [~,index_x_cut]=min(abs(cut_coord-x));
      x_cut=x(index_x_cut);
      fprintf('Extract variables at x=%.3f\n',x_cut);

    elseif strcmp(cut_plane,'y')

      [~,index_y_cut]=min(abs(cut_coord-y));
      y_cut=y(index_y_cut);
      fprintf('Extract variables at y=%.3f\n',y_cut);

    end
end

ntime=0;
for f = fileStart:fileStep:fileEnd
    ntime=ntime+1;
    fprintf('Timestep: %i \n',f);
    fname = sprintf('../../../restart/ruvwe.%07d.bin',f);
    fileID = fopen(fname);
    data = fread(fileID,ntot*6,'double');
    fclose(fileID);

    data=data(6:end);

    rho=data(1:ntot);
    u=data(ntot+1:2*ntot);
    v=data(2*ntot+1:3*ntot);
    w=data(3*ntot+1:4*ntot);
    e=data(4*ntot+1:5*ntot);

    rho = reshape(rho,imax,jmax,kmax);
    u = reshape(u,imax,jmax,kmax);
    v = reshape(v,imax,jmax,kmax);
    w = reshape(w,imax,jmax,kmax);
    e = reshape(e,imax,jmax,kmax);
   
if strcmp(visua,'2D') 
    if strcmp(cut_plane,'y')
        Rho_plane(ntime,:,:)=rho(:,index_cut,:);
        E_plane(ntime,:,:)=e(:,index_cut,:);
        W_plane(ntime,:,:)=w(:,index_cut,:);
        U_plane(ntime,:,:)=u(:,index_cut,:);
    elseif strcmp(cut_plane,'x')
        Rho_plane(ntime,:,:)=rho(index_cut,:,:);
        E_plane(ntime,:,:)=e(index_cut,:,:);
        W_plane(ntime,:,:)=w(index_cut,:,:);
        U_plane(ntime,:,:)=u(index_cut,:,:);
    elseif strcmp(cut_plane,'z')
        Rho_plane(ntime,:,:)=rho(:,:,index_cut);
        E_plane(ntime,:,:)=e(:,:,index_cut);
        W_plane(ntime,:,:)=w(:,:,index_cut);
        U_plane(ntime,:,:)=u(:,:,index_cut);
    end
    
elseif strcmp(visua,'1D')
    
       if strcmp(cut_plane,'z')

        Rho_cut(ntime,:,:)=rho(:,:,index_z_cut);
        E_cut(ntime,:,:)=e(:,:,index_z_cut);
        W_cut(ntime,:,:)=w(:,:,index_z_cut);
        U_cut(ntime,:,:)=u(:,:,index_z_cut);

    elseif strcmp(cut_plane,'x')

        Rho_cut(ntime,:,:)=rho(index_x_cut,:,:);
        E_cut(ntime,:,:)=e(index_x_cut,:,:);
        W_cut(ntime,:,:)=w(index_x_cut,:,:);
        U_cut(ntime,:,:)=u(index_x_cut,:,:);

    elseif strcmp(cut_plane,'y')

        Rho_cut(ntime,:,:)=rho(:,index_y_cut,:);
        E_cut(ntime,:,:)=e(:,index_y_cut,:);
        W_cut(ntime,:,:)=w(:,index_y_cut,:);
        U_cut(ntime,:,:)=u(:,index_y_cut,:);

       end
        
end

end
    %% Plot

if strcmp(visua,'2D')

    Rho_min=min(min(min(Rho_plane)));
    Rho_max=max(max(max(Rho_plane)));
    E_min=min(min(min(E_plane)));
    E_max=max(max(max(E_plane)));
    W_min=min(min(min(W_plane)));
    W_max=max(max(max(W_plane)));
    U_min=min(min(min(U_plane)));
    U_max=max(max(max(U_plane)));

    nsteps=200;

    if strcmp(cut_plane,'y')
        plot_coord1=Z;
        plot_coord2=X;
        coord2_Lim_max=10;
    elseif strcmp(cut_plane,'x')
        plot_coord1=Z;
        plot_coord2=Y;
        coord2_Lim_max=1;
    elseif strcmp(cut_plane,'z')
        plot_coord1=Y;
        plot_coord2=X;
        coord2_Lim_max=10;
    end

elseif strcmp(visua,'1D')

        if strcmp(cut_plane,'z')
            plot_coord=x;
            coord_max=5;
            index=jmax;
        elseif strcmp(cut_plane,'x')
            plot_coord=z;
            coord_max=z_max;
            index=jmax;
        elseif strcmp(cut_plane,'y')
            plot_coord=z;
            coord_max=z_max;
            index=imax;
        end
end

for i = 1:ntime  
    
    if strcmp(visua,'2D')

        Rho_plot(:,:)=Rho_plane(i,:,:);
        E_plot(:,:)=E_plane(i,:,:);
        W_plot(:,:)=W_plane(i,:,:);
        U_plot(:,:)=U_plane(i,:,:);
        
        figure(1)
        set(0, 'DefaultFigureRenderer', 'painters'); 
        subplot(2,2,1)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(plot_coord1,plot_coord2,Rho_plot,'LevelList',linspace(Rho_min-0.1,Rho_max+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        %xlabel(["$z/\delta_{99,in}$"])
        %ylabel("$x/\delta_{99,in}$")    
        ylim([0 coord2_Lim_max]);
        colormap(jet)
        cb.Label.String = '$\rho^*/\rho^*_{\infty}$';
        cb.Label.Interpreter = 'latex';

        subplot(2,2,2)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(plot_coord1,plot_coord2,E_plot,'LevelList',linspace(E_min-0.1,E_max+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        %xlabel("$x/\delta_{99,in}$")
        %ylabel("$y/\delta_{99,in}$")
        ylim([0 coord2_Lim_max]);
        colormap(jet)
        cb.Label.String = '$E^*/E^*_{\infty}$';
        cb.Label.Interpreter = 'latex';

        subplot(2,2,3)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(plot_coord1,plot_coord2,W_plot,'LevelList',linspace(0,W_max+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        %xlabel("$x/\delta_{99,in}$")
        %ylabel("$y/\delta_{99,in}$")
        ylim([0 coord2_Lim_max]);
        colormap(jet)
        cb.Label.String = '$u^*/u^*_{\infty}$';
        cb.Label.Interpreter = 'latex';

        subplot(2,2,4)
        c = parula(nsteps/2);
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex');
        colormap(c);  
        hold on
        [~,h] = contourf(plot_coord1,plot_coord2,U_plot*1e3,'LevelList',linspace(0,U_max*1e3+0.1,nsteps)); hold off
        set(h,'LineColor','none')
        ax=gca;
        ax.TickLabelInterpreter='latex';
        set(gca,'XMinorTick','on','YMinorTick','on')
        %xlabel("$x/\delta_{99,in}$")
        %ylabel("$y/\delta_{99,in}$")
        ylim([0 coord2_Lim_max]);
        colormap(jet)
        cb.Label.String = '$10^{3} \times v^*/u^*_{\infty} $';
        cb.Label.Interpreter = 'latex';
    
    elseif strcmp(visua,'1D')
         
        subplot(2,2,1)
        for i=1:index
            plot(W_cut(:,i),plot_coord,'-','LineWidth',3) 
        hold on
        end
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ u^*/u^*_\infty $','interpreter','latex');
        %ylabel("$y/\delta_{99,in}$");
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');
        axis([0 1.05 0 coord_max])

        subplot(2,2,2)
        for i=1:index
            plot(E_cut(:,i),plot_coord,'-','LineWidth',3) 
        hold on
        end
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ e^*/e^*_\infty $','interpreter','latex');
        axis([min(min(E_cut))-1 max(max(E_cut))+1 0 coord_max])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

        subplot(2,2,3)
        for i=1:index
            plot(Rho_cut(:,i),plot_coord,'-','LineWidth',3) 
        hold on
        end
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ \rho^*/\rho^*_\infty $','interpreter','latex');
        %ylabel("$y/\delta_{99,in}$");
        axis([0 1.05 0 coord_max])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');
        
        subplot(2,2,4)
        for i=1:index
            plot(U_cut(:,i),plot_coord,'-','LineWidth',3) 
        hold on
        end
        hold on
        grid off
        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlabel('$ 10^3 \times v^*/u^*_\infty $','interpreter','latex');
        ylim([0 coord_max])
        set(gca,'XMinorTick','on'); set(gca,'YMinorTick','on');

  
    end

end
