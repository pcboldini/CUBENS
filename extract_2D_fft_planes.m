%%%%%%%%%%%%%
clear all
clc
close all

mkdir postproc/planes_fft
path="postproc/planes";
path_new="postproc/planes_fft";

istart = 0000000;
iend = 0000100;
deltastep_old = 10;
deltastep_new = 10;
prop = ["r","u","v","w","e","p","t"]; % ,"mu","ka"
index_y=1;

%% Calculation

disp("Extracing planes:")
fprintf('steps: %d',istart);
fprintf(' - %d \n',iend);
fprintf('OLD delta_step: %d \n',deltastep_old);
fprintf('NEW delta_step: %d \n',deltastep_new);
fprintf('variable: %s \n',prop);

index=0;
for i = istart:deltastep_new:iend
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
    
    for j=1:numel(prop)
        fnames(index,j) = path +  "/ypl." + index_y + "." + prop(j) + "." + b + ".bin";
    end
end

index=0;
for i = istart:deltastep_old:iend
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
    y2(index) = b;
    
end

comp=contains(y,y2);
index=0;
for i = istart:deltastep_new:iend
   index=index+1;
   if islogical(comp(index))
      for j=1:numel(prop)
        file=fnames(index,j);
        disp(file)
        copyfile(file,path_new)
      end
   end
end

copyfile(path + "/x.bin",path_new);
copyfile(path + "/y.bin",path_new);
copyfile(path + "/z.bin",path_new);
disp(strcat("Extraction completed in ",path_new,"!"))
