%%%%%%%%%%%%%
clear all
clc
close all

istart = 0;
iend = 1000;
deltastep = 100;
var = 'w';
index_x=20;
index_y=1;
index_z=300;

show_plot='true'; % true/false

%%
disp("Extracting restart:")
fprintf('steps: %d',istart);
fprintf(' - %d \n',iend);
fprintf('Delta_step: %d \n',deltastep);
time_vec=(istart:deltastep:iend);

path="restart";
array_var=["rho","u","v","w","e"];
for i=1:length(array_var)
    for j=1:length(var)
        if array_var(i)==var(j)
            index_var(j)=i;
        end
    end
end

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
end

index=0;
for i = 1:numel(time_vec)
    index=index+1;

    fileID = fopen(fnames(index));
    data = fread(fileID,inf,'double');
    fclose(fileID);
    fprintf('Extract and reduce restart at %i \n',time_vec(i));

    imax=data(2);
    jmax=data(3);
    kmax=data(4);
    ntot=imax*jmax*kmax;

    data_head=data(1:5);
    data=data(6:end);
    
    for j=1:length(var)
        data_var=data(((index_var(j)-1)*ntot+1):(index_var(j)*ntot));
    end
  
    data_var_reshape=reshape(data_var,imax,jmax,kmax);
    data_plot(i)=data_var_reshape(index_x,index_y,index_z);

end

%% Plot

if strcmp(show_plot,'true')

    plot(time_vec,data_plot)
    xlabel("time")
    ylabel("amplitude")

end

