%%%%%%%%%%%%%
clear all
clc
close all

mkdir postproc/results/Yavg_select

path="postproc/results/Yavg";
path2="postproc";
path3="postproc/results";
path_planes="postproc/planes";
path_new="postproc/results/Yavg_select";

copyfile(path_planes + "/x.bin",path_new);
copyfile(path_planes + "/y.bin",path_new);
copyfile(path_planes + "/z.bin",path_new);

copyfile(path + "/Yave_r.bin",path_new);
copyfile(path + "/Yave_p.bin",path_new);
copyfile(path + "/Yave_T.bin",path_new);
copyfile(path + "/Yave_u.bin",path_new);
copyfile(path + "/Yave_v.bin",path_new);
copyfile(path + "/Yave_w.bin",path_new);
copyfile(path + "/Yave_mu.bin",path_new);
copyfile(path + "/Yave_T13.bin",path_new);
copyfile(path + "/Yave_q1.bin",path_new);
copyfile(path2 + "/wall_prop.txt",path_new);

copyfile(path3 + "/Re_tau.bin",path_new);
copyfile(path3 + "/u_tau.bin",path_new);
copyfile(path3 + "/Re_tau_sl.bin",path_new);
copyfile(path3 + "/u_tau_sl.bin",path_new);

disp(strcat("Extraction completed in ",path_new,"!"))

exit;
