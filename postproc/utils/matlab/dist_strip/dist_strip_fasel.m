clear all
clc
close all

set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Parameters

distz='classical'; % classical/alpha
dim='3D'; % 2D/3D
plots='T';
video='F';

z1=100.8291; % z_start
z2=200; % z_end
nz=1000;
ny=100;
Ly=0.6;
lamda_y=0.15;

fft_samples=50; % per period
fft_step=400; % per sample
n_periods=2; % number of periods
istart = 0;

amp1=0.5e-5;
amp2=5e-4;
omega1=0.0818;
omega2=0.10;
c_phase=0.298;

%% Calculation

alpha=omega1/c_phase;
lamda_z=2*pi/alpha;

Kmult=fft_samples*fft_step;
dt = 2*pi/omega1/Kmult;
freq=dt/2/pi;
samp_freq = (freq*fft_step)^-1; %Hz

iend=fft_samples*n_periods*fft_step;
ival= iend/fft_step+1;
freq_omega = samp_freq*(0:ival/2)/(ival-1);

t=linspace(0,dt*iend,ival);
y=linspace(0,Ly,ny);

if strcmp(distz,'classical')
    z=linspace(z1,z2,nz);
    zm=(z2+z1)/2;
    sigma=10;
elseif strcmp(distz,'alpha')
    z2=z1+lamda_z;
    z=linspace(z1,z2,nz);
    dz=z(2)-z(1);
    zm=(z2+z1)/2;
    sigma=0.03*lamda_z;

    freq_z=dz/2/pi;
    samp_z = (freq_z)^-1; 
    freq_alpha = samp_z*(0:nz/2)/(nz-1);
end

for i=1:numel(z)
    if z(i)<zm
        kappa=1;
        ksi(i)=(z(i)-z1)./(zm-z1);
    elseif z(i)>zm
        kappa=-1;
        ksi(i)=(z2-z(i))./(z2-zm);
    end
    
    f_z(i) = kappa*((15.1875*ksi(i).^5) - (35.4375*ksi(i).^4) + (20.25*ksi(i).^3));
    gaussian(i) = exp(-((z(i)-zm)./(2*sigma)).^2);
    
    ksi(i)=(z(i)-z1)./(zm-z1);
    gaussian_z(i) = exp(-((z(i)-zm)./(sigma)).^2);
    f_alpha_z(i) = sin(alpha*(z(i)-z1))*gaussian_z(i);
    ksi2(i)=(z(i)-zm)./(zm-z1);
    gaussian_z2(i) = exp(-(ksi2(i)./sigma).^2);
    f_alpha_z2(i) = sin(pi*ksi(i))*gaussian_z2(i);
end

if strcmp(dim,'2D')

for j=1:numel(t)
    u_wall_fk(j,:)=amp1*f_z*sin(omega1*t(j));
    u_wall_fl(j,:)=amp1*gaussian*sin(omega1*t(j));
    u_wall_alpha(j,:)=amp1*f_alpha_z2*sin(omega1*t(j));
    u_wall_test(j,:)=sin(omega1*t(j));
end

elseif strcmp(dim,'3D')

    for j=1:numel(y)
        g_y(j)=cos(2.0*pi*y(j)/lamda_y);
    end

    for k=1:numel(t)-1
        for j=1:numel(y)
            for i=1:numel(z)
                 u_wall1(k,j,i)=amp1*f_z(i)*sin(omega1*t(k));
                 u_wall2(k,j,i)=amp2*f_z(i)*g_y(j)*sin(omega2*t(k));
                 u_wall(k,j,i)=u_wall1(k,j,i)+u_wall2(k,j,i);
            end
        end
    end

end

%% FFT

% Y_fl = fft(u_wall_fl(1:end,:),[],1);
% P2_fl = abs(Y_fl/(ival));
% P1_fl = P2_fl(1:ival/2+1,:);
% P1_fl(2:end-1,:) = 2*P1_fl(2:end-1,:);
% 
% Y_fk = fft(u_wall_fk(1:end,:),[],1);
% P2_fk = abs(Y_fk/(ival));
% P1_fk = P2_fk(1:ival/2+1,:);
% P1_fk(2:end-1,:) = 2*P1_fk(2:end-1,:);
% 
% Y_test = fft(u_wall_test(1:end,:),[],1);
% P2_test = abs(Y_test/(ival));
% P1_test = P2_test(1:ival/2+1,:);
% P1_test(2:end-1,:) = 2*P1_test(2:end-1,:);
% 
% if strcmp(distz,'alpha')
%     Y_alpha_t = fft(u_wall_alpha(1:end,:),[],1);
%     P2_alpha_t = abs(Y_alpha_t/(ival));
%     P1_alpha_t = P2_alpha_t(1:ival/2+1,:);
%     P1_alpha_t(2:end-1,:) = 2*P1_alpha_t(2:end-1,:);
% 
%     Y_alpha_z = fft(u_wall_alpha(:,1:end),[],2);
%     P2_alpha_z = abs(Y_alpha_z/(nz));
%     P1_alpha_z = P2_alpha_z(:,1:nz/2+1);
%     P1_alpha_z(:,2:end-1) = 2*P1_alpha_z(:,2:end-1);
% end

%% Plots

if (plots=='T') 
    if strcmp(dim,'2D')

            for j=1:numel(t)
                figure(1)
                hold on
                plot(z,u_wall_alpha(j,:));
                hold off
                ylim([-amp1 amp1])
                pause(0.01)
            end

%             for i=1:numel(z)
%                 figure(2)
%                 %hold on
%                 plot(f,P1_fl(:,i));
%                 hold off
%             end

%             for j=1:numel(t)
%                 figure(2)
%                 hold on
%                 plot(freq_alpha,P1_alpha_z(j,:));
%                 hold off
%             end
    
        if (video=='T')
            movievector(t) = getframe(gcf);
        end
    
    elseif strcmp(dim,'3D')
    
        [Z,Y]=meshgrid(z,y);
    
      for k=1:numel(t)-1
    
            fig1=figure(1);
            %hold on
            %surf(Z,Y,squeeze(u_wall(k,:,:)));
            plot(z,squeeze(u_wall(k,1,:)));
            %plot(z,squeeze(u_wall(k,1,:)));
            hold off
            ylim([-amp2 amp2])
            set(findall(gcf,'type','text'),'Fontsize',26,'fontWeight','normal');
    
        if (video=='T')
            movievector(k) = getframe(gcf);
        end
    
      end

    end

end

%% Create video

if (video == 'T')

    MyWriter = VideoWriter('pert_p.mp4','MPEG-4');
    MyWriter.FrameRate = 10;
    open(MyWriter);
    writeVideo(MyWriter, movievector);
    close(MyWriter);  
end