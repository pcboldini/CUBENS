clear all
clc
close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Parameters

z1=41; % z_start
z2=51; % z_end

y1=0;
y2=9.63;

nz=1000;
ny=100;

amp1=2.25e-3;
amp2=7.5e-5;
omega1=0.193;
omega2=0.0963;
beta=1;

fft_samples=50; % per period
fft_step=400; % per sample
n_periods=2; % number of periods
istart = 0;

Kmult=fft_samples*fft_step;
dt = 2*pi/omega1/Kmult;
freq=dt/2/pi;
samp_freq = (freq*fft_step)^-1; %Hz

iend=fft_samples*n_periods*fft_step;
ival= iend/fft_step+1;
freq_omega = samp_freq*(0:ival/2)/(ival-1);

t=linspace(0,dt*iend,ival);

z=linspace(z1,z2,nz);
zm=(z2+z1)/2;
sigmaz=0.7;

y=linspace(y1,y2,ny);
ym=(y2+y1)/2;
sigmay=0.7;

for i=1:numel(z)
    gaussianz(i) = exp(-((z(i)-zm)./(2*sigmaz)).^2);
end

for i=1:numel(y)
    gaussiany(i) = 1 + 0.1*( exp(-((y(i)-ym-sigmay)./(sigmay)).^2)-exp(-((y(i)-ym+sigmay)./(sigmay)).^2) );
end

figure(1)
plot(z,gaussianz)
figure(2)
plot(y,gaussiany)


for k=1:numel(t)-1
    for j=1:numel(y)
        for i=1:numel(z)
             u_wall1(k,j,i)=gaussianz(i)*gaussiany(j);
             u_wall2(k,j,i)=amp1*sin(omega1*t(k))+amp2*sin(omega2*t(k)-beta*z(i));
             u_wall(k,j,i)=u_wall1(k,j,i)*u_wall2(k,j,i);
        end
    end
end

      for k=1:numel(t)-1
    
            fig1=figure(3);
            %hold on
            %surf(Z,Y,squeeze(u_wall(k,:,:)));
            plot(z,squeeze(u_wall(k,1,:)));
            %plot(z,squeeze(u_wall(k,1,:)));
            hold off
            ylim([-amp1 amp1])
            set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');
      end
