close all
clear all
clc

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontSize',26);

nz=2352;

lz=776.38;
xi=linspace(0,lz,nz)./lz;
dxi=xi(2)-xi(1);

fd1_coeff4=[1,-8,0,8,-1]./(12);
fd1_coeff1=[-1,1];

fd2_coeff4=[-1,16,-30,16,-1]./(12);
fd2_coeff2=[2,-5,4,-1];
fd2_coeff1=[1,-2,1];

% Transformation
% z=xi;
% zp=ones(nz,1)';
% zpp=zeros(nz,1)';
% dxi=lz*dxi;

zplus_min=10;
zplus_max=20;
z_1=0.25;
bump1=0.2;
delta1=z_1*bump1;
z_2=0.9;
bump2=0.04;
delta2=z_2*bump2;

z_int1=delta1*0.5*(zplus_min-zplus_max)*log(cosh((z_1)/delta1))...
      +delta2*0.5*(zplus_max-zplus_min)*log(cosh((z_2)/delta2));

for i=1:nz
    fact(i)=(i-1.0)/(nz-1.0);
    z(i)=(delta1*0.5*(zplus_min-zplus_max)*log(cosh((z_1-fact(i))/delta1))...
         +delta2*0.5*(zplus_max-zplus_min)*log(cosh((z_2-fact(i))/delta2))...
         +fact(i)*zplus_max-z_int1);
    zp(i)=(zplus_min+0.5*(zplus_max-zplus_min)*( 2-tanh((fact(i)-z_1)/(delta1))+tanh((fact(i)-z_2)/(delta2)) ) );
    zpp(i)=0.5*(zplus_max-zplus_min)*( (sech((z_2-fact(i))/delta2)).^2/delta2 - (sech((z_1-fact(i))/delta1)).^2/delta1 );
end

z_end=z(nz);
scaling=z_end/(lz);

z=z/scaling;
zp=zp/scaling;
zpp=zpp/scaling;

% Exact
y=sin(z/z(end)*pi);
dy=cos(z/z(end)*pi)/z(end)*pi;
ddy=-sin(z/z(end)*pi)/z(end)^2*pi^2;

% y=(z).^2;
% dy=2*(z);
% ddy=2.*ones(nz,1)';

% Approximation
dydz=zeros(nz-1,1);
dydxi=zeros(nz-1,1);

d2ydz2=zeros(nz-1,1);
d2ydxi2=zeros(nz-1,1);

zp_fd=zeros(nz,1);
zp_fd(1)=(z(2))/dxi;
for i=2:nz-1
    zp_fd(i)=(z(i+1)-z(i-1))/2/dxi;
end
zp_fd(end)=(z(end)-z(end-1))/dxi;

%% 1st derivative
% 4th-order (1st-forward at the boundaries)

for i=1:2
    for j=1:numel(fd1_coeff1)
        dydxi(i)=dydxi(i)+fd1_coeff1(j)*y(i+(j-1))/dxi;
    end
    dydz(i)=dydxi(i)./zp(i);
end

for i=3:nz-3
    for j=1:numel(fd1_coeff4)
        dydxi(i)=dydxi(i)+fd1_coeff4(j)*y(i+(j-3))/dxi;
    end
    dydz(i)=dydxi(i)./zp_fd(i);
end

for i=nz-2:nz-1
    for j=1:numel(fd1_coeff1)
        dydxi(i)=dydxi(i)+fd1_coeff1(j)*y(i+(j-1))/dxi;
    end
    dydz(i)=dydxi(i)./zp(i);
end

%% 2nd derivative

for i=3:nz-3
    for j=1:numel(fd2_coeff4)
        d2ydxi2(i)=d2ydxi2(i)+fd2_coeff4(j)*y(i+(j-3))/dxi^2;
    end
    d2ydz2(i)=1./(zp(i)).^2*d2ydxi2(i)-1./(zp(i)).^3*zpp(i)*dydxi(i);
end

for i=1:3
    for j=1:numel(fd2_coeff2)
        d2ydxi2(i)=d2ydxi2(i)+fd2_coeff2(j)*y(i+(j-1))/dxi^2;
    end
    if i>2
    else
        d2ydz2(i)=1./(zp(i)).^2*d2ydxi2(i)-1./(zp(i)).^3*zpp(i)*dydxi(i);
    end    
end

for i=nz-2:nz-1
    for j=1:numel(fd2_coeff1)
        d2ydxi2(i)=d2ydxi2(i)+fd2_coeff1(j)*y(i+(j-2))/dxi^2;
    end
    d2ydz2(i)=1./(zp(i)).^2*d2ydxi2(i)-1./(zp(i)).^3*zpp(i)*dydxi(i);
end

figure(1)
subplot(1,2,1)
plot(z,dy,z(1:end-1),dydz);
subplot(1,2,2)
plot(z,ddy,z(1:end-1),d2ydz2);