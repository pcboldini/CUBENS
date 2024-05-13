%%%%%%%%%%%%%
clear all
clc
close all

mkdir restart_new

disp("Extracing planes:")
istart = 100;
fprintf('steps: %d',istart);
iend = 100;
fprintf(' - %d \n',iend);
deltastep = 1;
fprintf('Delta_step: %d \n',deltastep);

path="restart";
path_planes="postproc/planes";
path_new="restart_new";

index=0;
for i = istart:deltastep:iend
    index=index+1;
    b = int2str(i);
    c = strlength(b);
    if istart>=1e6
        for k = 1:8-c
            b = "" + b;
        end
    else
        for k = 1:7-c
            b = "0" + b;
        end
    end
    y(index) = b;
    
    fnames(index) = path +  "/ruvwe." + b + ".bin";
    fnames_new(index) = path_new +  "/ruvwe." + b + ".bin";
end

index=0;
for i = istart:deltastep:iend
    index=index+1;
    b = int2str(i);

    fileID = fopen(fnames(index));
    data = fread(fileID,inf,'double');
    fclose(fileID);
    fprintf('Extract and reduce restart at %s \n',b);

    imax=data(2);
    jmax=data(3);
    kmax=data(4);
    ntot=imax*jmax*kmax;

    data_head=data(1:5);
    data=data(6:end);

    rho=data(1:ntot);
    u=data(ntot+1:2*ntot);
    w=data(3*ntot+1:4*ntot);

    interval=ntot;

    data_new(1:interval,1)=rho;
    data_new(interval+1:2*ntot)=u;
    data_new(2*interval+1:3*interval)=w;

    data_new=[data_head;data_new];

    fileID = fopen(fnames_new(index),'w');
    fwrite(fileID,data_new,'double');
    fclose(fileID);

end

copyfile(path_planes + "/x.bin",path_new);
copyfile(path_planes + "/y.bin",path_new);
copyfile(path_planes + "/z.bin",path_new);
disp(strcat("Extraction/reduction completed in ",path_new,"!"))
