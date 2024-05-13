%%%%%%%%%%%%%
clear all
clc
close all

%% Input

step=0000100;
jmax_new=10;

%% Proc

fname = sprintf('restart/ruvwe.%07d.bin',step);
fname2 = sprintf('restart/ruvwe.%07d_old.bin',step);
fileID = fopen(fname);
data = fread(fileID,inf,'double');
fclose(fileID);
fprintf('Extract restart at %i \n',step);

imax=data(2);
jmax=data(3);
kmax=data(4);
ntot=imax*jmax*kmax;
ntot_new=imax*jmax_new*kmax;

data_head=data(1:5);
data=data(6:end);

rho=data(1:ntot);
u=data(ntot+1:2*ntot);
v=data(2*ntot+1:3*ntot);
w=data(3*ntot+1:4*ntot);
e=data(4*ntot+1:5*ntot);

rho = reshape(rho,imax,jmax,kmax);
u = reshape(u,imax,jmax,kmax);
v = reshape(v,imax,jmax,kmax);
w = reshape(w,imax,jmax,kmax);
e = reshape(e,imax,jmax,kmax);

for j=1:jmax_new
    rho_new(:,j,:)=rho;
    u_new(:,j,:)=u;
    v_new(:,j,:)=v;
    w_new(:,j,:)=w;
    e_new(:,j,:)=e;
end

rho_new=reshape(rho_new,ntot_new,1);
u_new=reshape(u_new,ntot_new,1);
v_new=reshape(v_new,ntot_new,1);
w_new=reshape(w_new,ntot_new,1);
e_new=reshape(e_new,ntot_new,1);

interval=ntot_new;

data_head(3)=jmax_new;
data_new(1:interval,1)=rho_new;
data_new(interval+1:2*ntot_new)=u_new;
data_new(2*interval+1:3*interval)=v_new;
data_new(3*interval+1:4*interval)=w_new;
data_new(4*interval+1:5*interval)=e_new;

data_new=[data_head;data_new];

%% Write 

fileID = fopen(fname2,'w');
fwrite(fileID,data,'double');
fclose(fileID);

fileID = fopen(fname,'w');
fwrite(fileID,data_new,'double');
fclose(fileID);

disp("done!");




