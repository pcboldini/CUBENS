close all
clear all
clc

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26)

%fname =;
formatSpec = '%f';
fileID = fopen('../../xgrid.txt');
header = textscan(fileID, '%f',2);
header=header{:};
data = textscan(fileID, '%f%f%f%f','HeaderLines', 1,'CollectOutput',1);
data = data{:};
fileID = fclose(fileID);

ReTau=header(1);
stretch=header(2);
index=data(:,1);
x=data(:,2);
dx=data(:,3);
d2x=data(:,4);

%% Postprocessing

% actual mesh
disp('DNS mesh: ');
fprintf ('ReTau = %.3f \n',ReTau);
fprintf ('xStretch = %.3f \n',stretch);
index_max=index(end);
x_max=x(end);
fprintf ('nx = %.3f \n',index_max);
fprintf ('x_max = %.3f \n\n',x_max);

for i=1:numel(x)-1
Deltay(i)=x(i+1)-x(i);
end

% new mesh
ReTau_new = 50;
nx=300;
len_x=20;
stretch_new=3;

disp('new mesh:');
fprintf ('ReTau = %.3f \n',ReTau_new);
fprintf ('xStretch = %.3f \n',stretch_new);
fprintf ('nx = %.3f \n',nx);
fprintf ('x_max = %.3f \n',len_x);

yplFact = 0.6/ReTau_new*(nx-1)/len_x;

for i=1:nx
      fact(i)   =  (i-1.0)/(nx-1.0);
      tang(i) = (1 + tanh(stretch_new*(fact(i)-1)/2) / (tanh(stretch_new*0.5)));
      x_new(i)   = (fact(i)*yplFact + (1 + tanh(stretch_new*(fact(i)-1)/2) / (tanh(stretch_new*0.5)))*(1-yplFact))*len_x;
      dx_new(i)  = yplFact + (1-yplFact)*0.5*stretch_new/tanh(0.5*stretch_new)/cosh(stretch_new*(fact(i)-1)/2)^2.0;
      d2x_new(i) = -0.5*(1-yplFact)*(stretch_new^2)*tanh(stretch_new*(fact(i)-1)/2)/tanh(0.5*stretch_new)/(cosh(stretch_new*(fact(i)-1)/2))^2.0;      
end

index_new=linspace(1,numel(x_new),numel(x_new));
for i=1:numel(x_new)-1
    Deltay_new(i)=x_new(i+1)-x_new(i);
end

%% Plot

set(0, 'DefaultFigureRenderer', 'painters');
figure(1);
box
subplot(1,4,1)
plot(x,'-','LineWidth',3,'Color','blue') 
hold on
plot(x_new,'-','LineWidth',3,'Color','red') 
grid off
axis([0 300 0 x_max])
xlabel('Index','interpreter','latex');
ylabel('$y$','interpreter','latex');

subplot(1,4,2)
plot(index(1:end-1),Deltay,'-','LineWidth',3,'Color','blue') 
hold on
plot(index_new(1:end-1),Deltay_new,'-','LineWidth',3,'Color','red') 
grid off
xlim([0 300])
xlabel('Index','interpreter','latex');
ylabel('$\Delta y$','interpreter','latex');

subplot(1,4,3)
plot(index,dx,'-','LineWidth',3,'Color','blue') 
hold on
plot(index,dx_new,'-','LineWidth',3,'Color','red') 
grid off
xlim([0 index_max])
xlabel('Index','interpreter','latex');
ylabel('$\partial y$','interpreter','latex');

subplot(1,4,4)
plot(index,d2x,'-','LineWidth',3,'Color','blue') 
hold on
plot(index,d2x_new,'-','LineWidth',3,'Color','red') 
grid off
xlim([0 index_max])
xlabel('Index','interpreter','latex');
ylabel('$\partial^2 y$','interpreter','latex');




