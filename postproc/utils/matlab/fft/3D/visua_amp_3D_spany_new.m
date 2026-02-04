clear all
clc
close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Loadingb
path_case='../../../';

initBL=load(strcat(path_case,'initBL_param.mat')); 
params_variation=strcat(path_case,'params_variation.txt');
[initDNS]=import_DNS(params_variation);

% Parameters
Re_deltaSta=250;
Re_deltaEnd=350;
z_samp = [1 500];
istart = 0000; %1500000
iend = 1000; %1540000
istep = 100;
ival = (iend-istart)/istep+1;
prop = 'w'; % any var or ekin
x_crit='value'; % max/value
x_crit_val=0.2599; % value
refine_max='false'; % true/false
show_plot='true'; % true/false
save_plot='false'; % true/false

fprintf('steps: %d',istart);
fprintf('- %d \t',iend-istep);
fprintf('number of frames: %d \n',ival-1);

delta_eta_BL=initBL.BL_PARA.delta_eta_BL;
DNS_Re_99Sta=Re_deltaSta*delta_eta_BL;
z_start = Re_deltaSta^2/DNS_Re_99Sta;
F=initDNS.PERT_F;
Re_pert=initDNS.PERT_REMID;

fileID = fopen(strcat(path_case,'planes/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_case,'planes/y.bin'));
y = fread(fileID,'double');
fclose(fileID);
y = y - y(1);

fileID = fopen(strcat(path_case,'planes/z.bin'));
z = fread(fileID,'double');
fclose(fileID);
z = z - z(1);
dz = z(2)-z(1);

Re_plot=sqrt((z+z_start)*DNS_Re_99Sta);

nx = length(x);
ny = length(y);
nz = length(z);

nsamp=z_samp(2)-z_samp(1)+1;
gr_factor = Re_plot(z_samp(1):z_samp(2))./DNS_Re_99Sta;
gr_raw=zeros(1,nsamp);

index_harm=[1;2;3;4;5]; % 0.5,1,1.5,...
index_span=[0;1;2]; % 0,1,2,...
size_specry=numel(index_span);

specy=zeros(ival-1,nx,size_specry,nz); % 

count=0;

for f = istart:istep:iend-istep

count=count+1;

[specy] = load_spectra_new(path_case,prop,f,count,nx,index_span,nz,specy);

end

%% FFT

specy=specy(:,:,:,z_samp(1):z_samp(2));

[data_prop_ampt] = fourier_specy_new(specy,ival,index_harm);

%% Wall-normal treatment

[data_prop_amp] = calc_wallnormal(x_crit,x_crit_val,refine_max,nsamp,index_harm,data_prop_ampt,x);

%% Plot

% (omega=1/2*i,beta=j-1,k) omega for time and beta for spanwise, k is the streamwise

amp_0_0=squeeze(data_prop_amp(1,1,:));
amp_0_1=squeeze(data_prop_amp(1,2,:));
%amp_0_2=squeeze(data_prop_amp(1,3,:));

amp_12_0=squeeze(data_prop_amp(2,1,:));
amp_1_0=squeeze(data_prop_amp(3,1,:));
amp_2_0=squeeze(data_prop_amp(5,1,:));
%amp_3_0=squeeze(data_prop_amp(7,1,:));

amp_12_1=squeeze(data_prop_amp(2,2,:));
amp_1_1=squeeze(data_prop_amp(3,2,:));

for k = 2:nsamp-1
    gr_raw(k)= (amp_1_0(k+1) - amp_1_0(k-1))/2/dz;
end

alpha_i=gr_factor'.*gr_raw./amp_1_0';


if strcmp(show_plot,'true')

        xLim_min=(Re_pert);
        xLim_max=Re_deltaEnd;
        yLim_max=0;
        
        figure(1)
        box
        plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_1_0),'-','LineWidth',2, 'DisplayName',"$(F1, 0)$");
        hold on
        %plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_12_0),'--','LineWidth',2, 'DisplayName',"$(F1/2, 0)$");
        %plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_2_0),'-.','LineWidth',2, 'DisplayName',"$(F2, 0)$");
        %plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_3_0),'-.','LineWidth',2, 'DisplayName',"$(F3, 0)$");
        plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_12_1),'-','LineWidth',2, 'DisplayName',"$(F1/2, \beta)$");
        %plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_1_1),'--','LineWidth',2, 'DisplayName',"$(F1, \beta)$");
        %plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_0_0),'-','LineWidth',2, 'DisplayName',"$(F0, 0)$");
        plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_0_1),'-','LineWidth',2, 'DisplayName',"$(F0, \beta)$");
        %plot(Re_plot(z_samp(1):z_samp(2)),log10(amp_0_2),'-','LineWidth',2, 'DisplayName',"$(F0, 2\beta)$");

        ax=gca;
        ax.TickLabelInterpreter='latex';
        xlim([xLim_min xLim_max]);
        legend('boxoff','Location','southeast')
        ylim([-6 yLim_max]);
        xlabel('$ Re_\delta$','interpreter','latex');
        ylabel('$log_{10}(\tilde{q})$','interpreter','latex');
        set(gca,'Fontsize',26,'fontWeight','normal');
        set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');

end

if strcmp(save_plot,'true')

    DNS_RESULTS.Re_plot=Re_plot;
    DNS_RESULTS.z_samp=z_samp;
    DNS_RESULTS.Re_pert=Re_pert;
    DNS_RESULTS.data_prop_amp=data_prop_amp;

    saveFile         = num2str([initBL.BL_PARA.M_inf initBL.BL_PARA.T_inf initBL.BL_PARA.Ec_inf initBL.BL_PARA.Pr_inf initBL.BL_PARA.Re_deltaSta initBL.BL_PARA.Re_deltaEnd F nz nx], 'M%g_T%g_Ec%g_Pr%g_Re0%g-%g_F%g_Nz%g_Nx%g.mat');
    saveFile         = ['./results/DNS_AMP_' saveFile];
    save(saveFile,'initBL','initDNS','DNS_RESULTS');
    disp(['DNS saved in ' saveFile]);

end
