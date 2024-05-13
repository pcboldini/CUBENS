clear all
clc
% close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',20)


fileID = fopen('../../planes/x.bin');
x = fread(fileID,'double');
fclose(fileID);

fileID = fopen('../../planes/y.bin');
y = fread(fileID,'double');
fclose(fileID);

fileID = fopen('../../planes/z.bin');
z = fread(fileID,'double');
fclose(fileID);

imax = 1;
jmax = length(y);
kmax = length(z);

var = 'vorty';

fileStep  = 100;
fileStart = 0000000;
fileEnd   = 0002000-fileStep;

fname = sprintf('../../planes/xpl.%s.%07d.bin', var,fileStart);
fileID = fopen(fname);
data = fread(fileID,imax*jmax*kmax,'double');
fclose(fileID);

avg = data*0;
n = 0
for f = fileStart:fileStep:fileEnd
    n = n+1;
    fname = sprintf('../../planes/xpl.%s.%07d.bin', var,f);
    
    fileID = fopen(fname);
    data = fread(fileID,imax*jmax*kmax,'double');
    fclose(fileID);
    avg = avg+data;
end
avg = avg/n;


for f = fileStart:fileStep:fileEnd

    fname = sprintf('../../planes/xpl.%s.%07d.bin', var,f)
    
    fileID = fopen(fname);
    data = fread(fileID,imax*jmax*kmax,'double');
    fclose(fileID);

%    data=data-avg;
    
    data = reshape(data,imax,jmax,kmax);
    
    [X,Y] = meshgrid(z,x);

    data_plot(:,:)=data(:,1,:);

    contourf(X,Y,reshape(data(:,1,:),imax,kmax)); hold off
%    plot(z,data_plot(20,:)); hold on
    
    xlabel("z")
    ylabel("x")
 %   ylim([-1e-7 1e-7]);
 %    ylim([0 2e-6]);
    colormap(jet)
    pause(0.5)
end