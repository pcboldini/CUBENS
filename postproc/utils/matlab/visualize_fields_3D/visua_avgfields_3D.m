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

imax = length(x);
jmax = length(y);
kmax = length(z);

var = 'qx';

fname = sprintf('../../../restart/Yave_%s.bin', var);
fileID = fopen(fname);
data = fread(fileID,imax*jmax*kmax,'double');
fclose(fileID);

data = reshape(data,imax,1,kmax);

[X,Y] = meshgrid(z,x);
[C,h] =contourf(X,Y,reshape(data(:,1,:),imax,kmax),20); hold off
set(h,'LineColor','none')
xlabel("z")
ylabel("x")
colormap(jet)
