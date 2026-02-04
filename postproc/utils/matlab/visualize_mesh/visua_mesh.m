close all
clear all
clc

global slice

slice='y'; % z: streamwise, x: wall-normal, y: spanwise

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26)

path_planes="../../planes3D";

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

%% Postprocessing

if slice=='x'
    
    [Z,Y] = meshgrid(z,y);
    [Z Y];
    [Z' Y'];
    
elseif slice=='y'
    
    [Z,X] = meshgrid(z,x);
    [Z X];
    [Z' X'];
    
elseif slice=='z'

    [Y,X] = meshgrid(y,x);
    [Y X];
    [Y' X'];

end

%% Plot

index_zmax=max(z);
index_xmax=1;

if slice=='x'
    
    set(0, 'DefaultFigureRenderer', 'painters');
    figure(1);
    box
    plot(Z,Y,'k');
    hold on
    plot(Z',Y','k');
%    xlim([0 index_zmax]);
%    ylim([0 index_xmax]);
    xlabel('$x/\delta_{99,start}$','interpreter','latex');
    ylabel('$z/\delta_{99,start}$','interpreter','latex');
    
elseif slice=='y'
    
    set(0, 'DefaultFigureRenderer', 'painters');
    figure(1);
    box
    plot(Z,X,'k');
    hold on
    plot(Z',X','k');
    xlim([0 index_zmax]);
    ylim([0 index_xmax]);
    xlabel('$x/\delta_{99,start}$','interpreter','latex');
    ylabel('$y/\delta_{99,start}$','interpreter','latex');
    
elseif slice=='z'
    
    set(0, 'DefaultFigureRenderer', 'painters');
    figure(1);
    box
    plot(Y,X,'k');
    hold on
    plot(Y',X','k');
%     xlim([0 index_zmax]);
%     ylim([0 index_xmax]);
    xlabel('$z/\delta_{99,start}$','interpreter','latex');
    ylabel('$y/\delta_{99,start}$','interpreter','latex');
    
end



