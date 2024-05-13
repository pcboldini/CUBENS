clear all
clc
close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Loading
path_case="../../../";
path_planes="";
SD='';
GR='';
initBL='postproc/initBL_param.mat'; 
init_params='postproc/params_variation.txt';

show_plot='true'; % true/false
save_plot='false'; % true/false

% Parameters
Re_deltaSta=250;
Re_deltaEnd=350;
z_samp = [20 2500];
istart = 4000000;
iend = 4060000;
istep = 600;
ival = 1+(iend-istart)/istep;
prop = 'w';
x_crit='max'; % max/value/ekin
x_crit_val=0.25; % value
nsmooth = 10;

fprintf('steps: %d',istart);
fprintf('- %d \t',iend-istep);
fprintf('number of frames: %d \n',ival-1);

%%

[initDNS]=import_DNS(init_params,path_case);
initBL=load(strcat('../',path_case,'/',initBL)); 
% SD=load(strcat('../',path_case,'/',SD));
% GR=load(strcat('../',path_case,'/',GR));

delta_eta_BL=initBL.BL_PARA.delta_eta_BL;
DNS_Re_99Sta=Re_deltaSta*delta_eta_BL;
z_start = Re_deltaSta^2/DNS_Re_99Sta;
F=initDNS.PERT_F;
Re_pert=initDNS.PERT_REMID;

fileID = fopen(strcat('../',path_case,'/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat('../',path_case,'/y.bin'));
y = fread(fileID,'double');
fclose(fileID);
y = y - y(1);

fileID = fopen(strcat('../',path_case,'/z.bin'));
z = fread(fileID,'double');
fclose(fileID);
z = z - z(1);
dz = z(2)-z(1);

Re_plot=sqrt((z+z_start)*DNS_Re_99Sta);

nx = length(x);
ny = 1;
nz = length(z);

nsamp=z_samp(2)-z_samp(1)+1;
gr_factor = Re_plot(z_samp(1):z_samp(2))./DNS_Re_99Sta;
cr_factor = DNS_Re_99Sta*F;
gr_raw=zeros(1,nsamp);
phase_raw=zeros(1,nsamp);

count=0;
for f = istart:istep:(iend-istep)
    count=count+1;

    if strcmp(x_crit,'ekin')

        fname=sprintf(strcat('../',path_case,'/',path_planes,'/ypl.%s.%07d.bin'), 'w',f);
        fileID = fopen(fname);
        data_w = fread(fileID,nx*ny*nz,'double');
        fclose(fileID);

        data_w = reshape(data_w,nx,1,nz);
        data_w_cut=data_w(:,:,z_samp(1):z_samp(2));
        data_w_time(count,:,:) = data_w_cut(:,1,:);
        
        fname=sprintf(strcat('../',path_case,'/',path_planes,'/ypl.%s.%07d.bin'), 'u',f);
        fileID = fopen(fname);
        data_u = fread(fileID,nx*ny*nz,'double');
        fclose(fileID);

        data_u = reshape(data_u,nx,1,nz);
        data_u_cut=data_u(:,:,z_samp(1):z_samp(2));
        data_u_time(count,:,:) = data_u_cut(:,1,:);

        fname=sprintf(strcat('../',path_case,'/',path_planes,'/ypl.%s.%07d.bin'), 'r',f);
        fileID = fopen(fname);
        data_r = fread(fileID,nx*ny*nz,'double');
        fclose(fileID);

        data_r = reshape(data_r,nx,1,nz);
        data_r_cut=data_r(:,:,z_samp(1):z_samp(2));
        data_r_time(count,:,:) = data_r_cut(:,1,:);

    else

        fname=sprintf(strcat('../',path_case,'/',path_planes,'/ypl.%s.%07d.bin'), prop,f);
        fileID = fopen(fname);
        data = fread(fileID,nx*ny*nz,'double');
        fclose(fileID);
    
        data = reshape(data,nx,1,nz);
        data_cut=data(:,:,z_samp(1):z_samp(2));
        data_time(count,:,:) = data_cut(:,1,:);

    end

    fname_p=sprintf(strcat('../',path_case,'/',path_planes,'/ypl.%s.%07d.bin'), 'p',f);
    fileID_p = fopen(fname_p);
    data_p = fread(fileID_p,nx*ny*nz,'double');
    fclose(fileID_p);

    data_p = reshape(data_p,nx,1,nz);
    data_p_cut=data_p(:,:,z_samp(1):z_samp(2));
    data_p_time(count,:,:) = data_p_cut(:,1,:);
end

if strcmp(x_crit,'ekin')

    data_w_mean = mean(data_w_time);
    data_w_fluc=data_w_time-data_w_mean;

    data_w_fft = fft(data_w_fluc,[],1);
    data_w_fft2 = abs(data_w_fft/(ival-1));

    data_w_fft_abs = data_w_fft2(1:round(ival/2,0),:,:);
    data_w_fft_abs(2:end-1,:,:) = 2*data_w_fft_abs(2:end-1,:,:);

    data_u_mean = mean(data_u_time);
    data_u_fluc=data_u_time-data_u_mean;

    data_u_fft = fft(data_u_fluc,[],1);
    data_u_fft2 = abs(data_u_fft/(ival-1));

    data_u_fft_abs = data_u_fft2(1:round(ival/2,0),:,:);
    data_u_fft_abs(2:end-1,:,:) = 2*data_u_fft_abs(2:end-1,:,:);

    data_r_mean = mean(data_r_time);

    [vals_w_time,idxs_w_time1] = max(data_w_fft_abs,[],1);
    [vals_u_time,idxs_u_time1] = max(data_u_fft_abs,[],1);

else

    data_mean = mean(data_time);
    data_fluc=data_time-data_mean;
    
    data_fft = fft(data_fluc,[],1);
    data_fft2 = abs(data_fft/(ival-1));
    
    data_fft_abs = data_fft2(1:round(ival/2,0),:,:);
    data_fft_abs(2:end-1,:,:) = 2*data_fft_abs(2:end-1,:,:);

    [vals_time,idxs_time1] = max(data_fft_abs,[],1);
    [~,idxs_x] = max(vals_time,[],2);
    idxs_time=idxs_time1(idxs_x);

end

data_p_mean = mean(data_p_time);
data_p_fluc=data_p_time-data_p_mean;
data_p_fft = fft(data_p_fluc,[],1);
data_phase2 = angle(data_p_fft);

data_phase = data_phase2(1:round(ival/2,0),:,:);
data_phase(2:end-1,:,:) = 2*data_phase(2:end-1,:,:);


if strcmp(x_crit,'max')
    for k=1:nsamp
        % Refinement around wall-normal max value
        xmax_fine=linspace(x(idxs_x(k)-3),x(idxs_x(k)+3),numel(x));
        vals_time_fine=spline(x,vals_time(1,:,k),xmax_fine);

        data_amp(k)=max(vals_time_fine);
        data_phase_max(k)=data_phase2(idxs_time(k),idxs_x(k),k);
    end
elseif strcmp(x_crit,'ekin')
        ekin(:,:)=0.5*data_r_mean.*(vals_w_time.^2+vals_u_time.^2);
        ekin_int(:)=sqrt(trapz(x,ekin(:,:),1));
        data_amp=ekin_int;

        [~,idxs_x] = max(ekin,[],1);
        idxs_time=idxs_w_time1(idxs_x);
        for k=1:nsamp
            data_phase_max(k)=data_phase2(idxs_time(k),idxs_x(k),k);
        end
elseif strcmp(x_crit,'value')
    [~,index_x_crit_val]=min(abs(x_crit_val-x)); 
    for k=1:nsamp
        data_amp(k)=max(data_fft_abs(:,index_x_crit_val,k));
        data_phase_max(k)=data_phase2(idxs_time(k),index_x_crit_val,k);
    end
end

 data_phase_max = unwrap(data_phase_max);

for k = 2:nsamp-1
    gr_raw(k)= (data_amp(k+1) - data_amp(k-1))/2/dz;
    phase_raw(k)= (data_phase_max(k+1) - data_phase_max(k-1))/2/dz;
end

alpha_i=gr_factor'.*gr_raw./data_amp;
alpha_r=gr_factor'.*abs(phase_raw);
cr=cr_factor'./phase_raw;

% for i = nsmooth:length(alpha_i)-nsmooth+1
%     alpha_i_smooth(i) = sum(alpha_i(i-nsmooth+1:i+nsmooth-1))/(2*nsmooth-1);
%     alpha_r_smooth(i) = sum(alpha_r(i-nsmooth+1:i+nsmooth-1))/(2*nsmooth-1);
%     cr_smooth(i) = sum(cr(i-nsmooth+1:i+nsmooth-1))/(2*nsmooth-1);
% end
% 
% plot_vec = (1:20:length(alpha_i_smooth));


%% Plot

if strcmp(show_plot,'true')

alphaVal_max=0.007;
alphaVal_num=30;

cVal_max=0.6;
cVal_num=10;

xLim_max=2000;
xLim_min=0;
yLim_max=300;

Fmin=min(SD.LST_PARA.F_all);
Fmax=max(SD.LST_PARA.F_all);    
noF=length(SD.LST_PARA.F_all); 
F_r_StabDiag=SD.LST_PARA.F_all;
noRe0=length(SD.LST_PARA.Re0_all);

Eval_all_imag=imag(SD.LST_RESULTS.Eval_all);
Eval_all_real=real(SD.LST_RESULTS.Eval_all);

F_StabDiag=linspace(Fmin, Fmax, 1*noF);
Re0_StabDiag=linspace(min(SD.LST_PARA.Re0_all),max(SD.LST_PARA.Re0_all),1*noRe0);

[Rex_axis F_axis] = meshgrid(Re0_StabDiag, F_StabDiag);
ALPHA_i_StabDiag = griddata(Re0_StabDiag,F_r_StabDiag,Eval_all_imag,Rex_axis,F_axis,'cubic');

SD.LST_PARA.omega_all=SD.LST_PARA.F_all.*SD.LST_PARA.Re0_all;
SD.LST_PARA.c_r=SD.LST_PARA.omega_all'./Eval_all_real;

GR.LST_PARA.omega_all=GR.LST_PARA.F_0.*GR.LST_PARA.Re0_all;
GR.LST_PARA.c_r=GR.LST_PARA.omega_all'./real(GR.LST_RESULTS.Eval_all);

c = jet(15);
figure(1)
subplot(2,2,1)
box
cb=colorbar;
set(cb,'TickLabelInterpreter','latex');
colormap(c);
hold on
contourf(Rex_axis,F_axis*1e6,-ALPHA_i_StabDiag,'LineWidth',1,'LineColor',[0 0 0],'LevelList',linspace(0,alphaVal_max,alphaVal_num));
line([Re_deltaSta Re_deltaEnd],[F*1e6 F*1e6],'LineStyle','--','LineWidth',2,'Color',[0 0 0]); % 
plot(Re_pert,F*1e6,'p','LineWidth',3,'Color','red','MarkerSize',10,'MarkerFaceColor','black')
ax=gca;
ax.CLim=[0 alphaVal_max];
ax.TickLabelInterpreter='latex';
xlim([xLim_min xLim_max]);
ylim([0 yLim_max]);
xlabel('$ Re_\delta$','interpreter','latex');
ylabel('$ F \times 10^6$','interpreter','latex');
set(gca,'Fontsize',24,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',24,'fontWeight','normal');

subplot(2,2,2)
plot(Re_plot(z_samp(1):z_samp(2)),alpha_i,'-','Color','blue')
hold on
plot(GR.LST_PARA.Re0_all,-imag(GR.LST_RESULTS.Eval_all),'LineWidth',2,'Color', [0 0 0])
yline(0)
ylim([-0.02 0.02])
xlim([600 1700])
legend('DNS','LST','Location','north')
ylabel('Growth Rate $\alpha_i$')
xlabel('$Re_\mathrm{\delta}$')
set(gca,'Color',[1 1 1],'Fontsize',20);
set(gca,'Fontsize',24,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',24,'fontWeight','normal');

subplot(2,2,3)
plot(Re_plot(z_samp(1):z_samp(2)),alpha_r,'-','Color','blue')
hold on%
plot(GR.LST_PARA.Re0_all,real(GR.LST_RESULTS.Eval_all),'LineWidth',1,'Color', [0 0 0])
yline(0)
%ylim([0.2 0.6])
xlim([600 1700])
legend('DNS','LST','Location','north')
ylabel('Streamwise wavenumber $\alpha_r$')
xlabel('$Re_\mathrm{\delta}$')
set(gca,'Color',[1 1 1],'Fontsize',20);
set(gca,'Fontsize',24,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',24,'fontWeight','normal');

subplot(2,2,4)
plot(Re_plot(z_samp(1):z_samp(2)),abs(cr),'-','Color','blue')
hold on
plot(GR.LST_PARA.Re0_all,GR.LST_PARA.c_r,'LineWidth',1,'Color', [0 0 0])
yline(0)
ylim([0.2 0.6])
xlim([600 1700])
legend('DNS','LST','Location','north')
ylabel('Phase speed $c_r$')
xlabel('$Re_\mathrm{\delta}$')
set(gca,'Color',[1 1 1],'Fontsize',20);
set(gca,'Fontsize',24,'fontWeight','normal');
set(findall(gcf,'type','text'),'Fontsize',24,'fontWeight','normal');

end

if strcmp(save_plot,'true')

    DNS_RESULTS.plot_vec=plot_vec;
    DNS_RESULTS.Re_plot=Re_plot;
    DNS_RESULTS.F=F;
    DNS_RESULTS.z_samp=z_samp;
    DNS_RESULTS.Re_pert=Re_pert;
    DNS_RESULTS.alpha_i=alpha_i;
    DNS_RESULTS.alpha_i_smooth=alpha_i_smooth;
    DNS_RESULTS.cr_smooth=cr_smooth;

    saveFile         = num2str([initBL.BL_PARA.M_inf initBL.BL_PARA.T_inf initBL.BL_PARA.Ec_inf initBL.BL_PARA.Pr_inf initBL.BL_PARA.Re_deltaSta initBL.BL_PARA.Re_deltaEnd F nz nx], 'M%g_T%g_Ec%g_Pr%g_Re0%g-%g_F%g_Nz%g_Nx%g.mat');
    saveFile         = ['./Results/DNS_GR_' saveFile];
    save(saveFile,'initBL','initDNS','GR','SD','DNS_RESULTS');
    disp(['DNS saved in ' saveFile]);

end


