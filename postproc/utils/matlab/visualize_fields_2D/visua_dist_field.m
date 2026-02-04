clear all
clc
% close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',20)

show_plot='cut'; % contour/cut
path_planes="../../planes";

fileID = fopen(strcat(path_planes,'/x.bin'));
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_planes,'/y.bin'));
y = fread(fileID,'double');
fclose(fileID);

fileID = fopen(strcat(path_planes,'/z.bin'));
z = fread(fileID,'double');
fclose(fileID);

imax = length(x);
jmax = length(y);
kmax = length(z);

var = 'r';
video = 'F';

fileStep  = 100;
fileStart = 0000000;
fileEnd   = 0000100; % fileStep

fname = sprintf(strcat(path_planes,'/ypl.1.%s.%07d.bin'), var,fileStart);
fileID = fopen(fname);
data = fread(fileID,imax*jmax*kmax,'double');
fclose(fileID);

avg = data*0;
n = 0;
for f = fileStart:fileStep:fileEnd
    n = n+1;
    fname = sprintf(strcat(path_planes,'/ypl.1.%s.%07d.bin'), var,f);
    
    fileID = fopen(fname);
    data = fread(fileID,imax*jmax*kmax,'double');
    fclose(fileID);
    avg = avg+data;
end
avg = avg/n;

n = 0;
for f = fileStart:fileStep:fileEnd
    n = n+1;
    fname = sprintf(strcat(path_planes,'/ypl.1.%s.%07d.bin'), var,f);
    
    fileID = fopen(fname);
    data = fread(fileID,imax*jmax*kmax,'double');
    fclose(fileID);

    data=data-avg;
    
    data2 = reshape(data,imax,jmax,kmax);
    
    [X,Y] = meshgrid(z,x);

    data_plot(:,:)=data2(:,1,:);

    max_data=max(max(data_plot));

    if strcmp(show_plot,'contour')

        fig1=figure(1);
        box
        colorbar('delete') 
        cb=colorbar;
        cb.Label.String = '$u''/u''_{max}$';
        cb.Location = 'northoutside';
        cb.Label.Interpreter = 'latex';
        set(cb,'TickLabelInterpreter','latex')
        colormap(jet)
        hold on
        contourf(X,Y,reshape(data(:,1,:),imax,kmax)./1); hold off  
        ax=gca;
        xlabel('$ z/\delta_{99,in}$','interpreter','latex');
        ylabel('$ x/\delta_{99,in}$','interpreter','latex');
        ax.TickLabelInterpreter='latex';
        set(gca,'Fontsize',26,'fontWeight','normal');
        set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');
        set(gcf, 'Position', get(0, 'Screensize'));
        pause(0.1)

    elseif strcmp(show_plot,'cut')

        fig1=figure(1);
        hold on
        plot(z,data_plot(20,:)); 
        hold off
        pause(0.1)

    end

    if (video=='T')
%         ax.Units = 'pixels';
%         pos = ax.Position;
%         ti = ax.TightInset;
%         rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        movievector(n) = getframe(gcf);
    end

end


%% Create video

if (video == 'T')

    MyWriter = VideoWriter('pert_p.mp4','MPEG-4');
    MyWriter.FrameRate = 50;
    open(MyWriter);
    writeVideo(MyWriter, movievector);
    close(MyWriter);  
end