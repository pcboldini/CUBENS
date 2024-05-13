%%%%%%%%%%%%%
clear all
clc
close all

mkdir postproc/results/fft_select

disp("Extracing fft:")
istart = 000000;
fprintf('steps: %d',istart);
iend = 0000010;
fprintf(' - %d \n',iend);
deltastep_old = 1;
fprintf('OLD delta_step: %d \n',deltastep_old);
deltastep_new = 1;
fprintf('NEW delta_step: %d \n',deltastep_new);

var={'r','u','w'};
index_fft_span = [0,1,2];

path="postproc/results/fft";
path_planes="postproc/planes";
path_new="postproc/results/fft_select";

index=0;
for i = istart:deltastep_new:(iend-1)
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
    
    for jj=1:numel(var)
        for j=1:numel(index_fft_span)
        fnames(index,j,jj) = path +  "/fft_" + var{jj} +"_span" + index_fft_span(j) + '_' + b + ".bin";
        end
    end
end

index=0;
for i = istart:deltastep_old:(iend-1)
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

     for jj=1:numel(var)
        for j=1:numel(index_fft_span)
            y2(index,j,jj) = b;
        end
     end
    
end

comp=contains(y,y2);
index=0;
for i = istart:deltastep_new:(iend-1)
   index=index+1;
   if islogical(comp(index))
       for jj=1:numel(var)
        for j=1:numel(index_fft_span)
            file=fnames(index,j,jj);
            disp(file)
            copyfile(file,path_new)
        end
       end
   end
end

copyfile(path_planes + "/x.bin",path_new);
copyfile(path_planes + "/y.bin",path_new);
copyfile(path_planes + "/z.bin",path_new);
disp(strcat("Extraction completed in ",path_new,"!"))

exit;

