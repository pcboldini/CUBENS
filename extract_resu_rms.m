%%%%%%%%%%%%%
clear all
clc
close all

mkdir postproc/results/rms_select

path="postproc/results/rms";
path2="postproc/results";
path_planes="postproc/planes";
path_new="postproc/results/rms_select";

copyfile(path_planes + "/x.bin",path_new);
copyfile(path_planes + "/y.bin",path_new);
copyfile(path_planes + "/z.bin",path_new);

copyfile(path2 + "/Re_tau.bin",path_new);
copyfile(path2 + "/u_tau.bin",path_new);
copyfile(path2 + "/Re_tau_sl.bin",path_new);
copyfile(path2 + "/u_tau_sl.bin",path_new);
copyfile(path + "/aProd_wfluc_wfluc.bin",path_new);
copyfile(path + "/aProd_wfluc_ufluc.bin",path_new);
copyfile(path + "/aProd_ufluc_ufluc.bin",path_new);
copyfile(path + "/aProd_vfluc_vfluc.bin",path_new);
copyfile(path + "/aProd_rfluc_rfluc.bin",path_new);
copyfile(path + "/aProd_mufluc_mufluc.bin",path_new);
copyfile(path + "/aProd_kafluc_kafluc.bin",path_new);
copyfile(path + "/aProd_Cpfluc_Cpfluc.bin",path_new);


disp(strcat("Extraction completed in ",path_new,"!"))

exit;
