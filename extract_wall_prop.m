%%%%%%%%%%%%%
clear all
clc
close all

mkdir postproc/wall_prop

path="postproc/results";
path_planes="postproc/planes";
path_new="postproc/wall_prop";

copyfile(path_planes + "/x.bin",path_new);
copyfile(path_planes + "/y.bin",path_new);
copyfile(path_planes + "/z.bin",path_new);

copyfile(path + "/Yave_T13.bin",path_new);
copyfile(path + "/Yave_q1.bin",path_new);

disp(strcat("Extraction completed in ",path_new,"!"))
