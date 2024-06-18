clear all
clc
close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26)

%% Parameters 

prop = 'U';
USE_EOS='VdW'; % IG/VdW
video = 'F';
path_new="../../planes";

istart = 0000000;
iend = 0000100;
deltastep = 10;

%% Calculation

fprintf('steps: %d',istart);
fprintf(' - %d \n',iend);
fprintf('OLD delta_step: %d \n',deltastep);

if strcmp(USE_EOS,'IG')
initBL=load('../../initBL_param_IG.mat');
elseif strcmp(USE_EOS,'VdW')
initBL=load('../../initBL_param.mat');
end

Re_deltaSta=initBL.BL_PARA.Re_deltaSta;
Re_deltaEnd=initBL.BL_PARA.Re_deltaEnd;
DNS_Re_99Sta=initBL.DNS_PARA.DNS_Re_99Sta;
z_start=initBL.DNS_PARA.DNS_x_Sta;

fileID = fopen(strcat(path_new,'/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_new,'/z.bin'));
z = fread(fileID,'double');
fclose(fileID);
z = z - z(1);

BL_scaling=sqrt((z(end)+initBL.DNS_PARA.DNS_x_Sta)./initBL.DNS_PARA.DNS_x_Sta);
diff_bl=BL_scaling-x; [~,index_BL_scaling]=min(abs(diff_bl));
x99=x(index_BL_scaling);

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
for i = istart:deltastep:iend
    ii = (i - istart)/deltastep + 1;
    b = int2str(i);
    c = strlength(b);
    for k = 1:7-c
        b = "0" + b;
    end
    fname = path_new + "/ypl." + prop + "." + b + ".bin";
    
    fileID = fopen(fname);
    mat_read = fread(fileID,[nx nz],'double');
    fclose(fileID);

    for k = 1:nz
    A_plot(:,k)=mat_read(:,k);
    end

    A_max=max(max(A_plot));
    A_min=min(min(A_plot));
   
    A_limit=[0.99*A_min 1.01*A_max];
    
    xlabel('$z/\delta_\mathrm{in}$')
    ylabel('$x/\delta_\mathrm{in}$')
    hold on
    contour(z(1:nz),x(1:nx),A_plot,linspace(A_min,A_max,30))
    shading interp
    colorbar
    hold off
    ylim([0 x99*2]);
    title(prop);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    pause(0.1)
    
    if (video=='T')
        F(ii) = getframe(gcf); 
    end

end

%% Create video

if (video == 'T')
    writerObj = VideoWriter('../../results/base_transient_2D.mp4','MPEG-4');
    writerObj.FrameRate = 5;

    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end
