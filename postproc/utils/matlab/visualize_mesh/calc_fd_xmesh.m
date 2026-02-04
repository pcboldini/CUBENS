close all
clear all
clc

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26);

nx=100000;

lx=20;
xi=linspace(0,lx,nx)/(lx);
dxi=xi(2)-xi(1);

fd1_coeff4=[1,-8,0,8,-1]./(12);
fd1_coeff1=[-1,1];

fd2_coeff4=[-1,16,-30,16,-1]./(12);
%fd2_coeff4=[-30]./(12);
fd2_coeff2=[2,-5,4,-1];
fd2_coeff1=[1,-2,1];

% Transformation

x=xi;
xp=ones(nx,1)';
xpp=zeros(nx,1)';
dxi=dxi;

ReTau=44.472;
yplFact = 1.3484769035831690;
stretch=5;

for i=1:nx
      fact   =  (i-1.0)/(nx-1.0);
      x(i)   = (fact*yplFact + (1.0 + tanh(stretch*(fact-1.0)/2.0) / (tanh(stretch*0.5)))...
                *(1.0-yplFact))*lx;
      xp(i)  = (yplFact + (1-yplFact)*0.5*stretch/tanh(0.5*stretch)/cosh(stretch*(fact-1)/2)^2)*lx;
      xpp(i) = (-0.5*(1-yplFact)*(stretch^2)*tanh(stretch*(fact-1)/2)/tanh(0.5*stretch)/(cosh(stretch*(fact-1)/2))^2.0)*lx;
end

% Exact
y=sin(x/x(end)*pi);
dy=cos(x/x(end)*pi)*pi/x(end);
ddy=-sin(x/x(end)*pi)*pi^2/x(end)/x(end);

% y=(x).^2;
% dy=2*(x)*1;
% ddy=2.*ones(nx,1)';

% Approximation
dydz=zeros(nx-1,1);
dydxi=zeros(nx-1,1);

d2ydz2=zeros(nx-1,1);
d2ydxi2=zeros(nx-1,1);

%% 1st derivative
% 4th-order (1st-forward at the boundaries)

for i=1:2
    for j=1:numel(fd1_coeff1)
        dydxi(i)=dydxi(i)+fd1_coeff1(j)*y(i+(j-1))/dxi;
    end
    dydz(i)=dydxi(i)./xp(i);
end

for i=3:nx-3
    for j=1:numel(fd1_coeff4)
        dydxi(i)=dydxi(i)+fd1_coeff4(j)*y(i+(j-3))/dxi;
    end
    dydz(i)=dydxi(i)./xp(i);
end

for i=nx-2:nx-1
    for j=1:numel(fd1_coeff1)
        dydxi(i)=dydxi(i)+fd1_coeff1(j)*y(i+(j-1))/dxi;
    end
    dydz(i)=dydxi(i)./xp(i);
end

%% 2nd derivative

for i=3:nx-3
    for j=1:numel(fd2_coeff4)
        d2ydxi2(i)=d2ydxi2(i)+fd2_coeff4(j)*y(i+(j-3))/dxi^2;
    end
    d2ydz2(i)=1./(xp(i)).^2*d2ydxi2(i)-1./(xp(i)).^3*xpp(i)*dydxi(i);
end

for i=1:3
    for j=1:numel(fd2_coeff2)
        d2ydxi2(i)=d2ydxi2(i)+fd2_coeff2(j)*y(i+(j-1))/dxi^2;
    end
    if i>2
    else
        d2ydz2(i)=1./(xp(i)).^2*d2ydxi2(i)-1./(xp(i)).^3*xpp(i)*dydxi(i);
    end    
end

for i=nx-2:nx-1
    for j=1:numel(fd2_coeff1)
        d2ydxi2(i)=d2ydxi2(i)+fd2_coeff1(j)*y(i+(j-2))/dxi^2;
    end
    d2ydz2(i)=1./(xp(i)).^2*d2ydxi2(i)-1./(xp(i)).^3*xpp(i)*dydxi(i);
end

figure(1)
subplot(1,2,1)
plot(x,dy,x(1:end-1),dydz);
subplot(1,2,2)
plot(x,ddy,x(1:end-1),d2ydz2);
%ylim([-1.9 2.1])