import shutil
import os
import numpy as np
from writexmf import writexmf

current_path=os.getcwd()
dst_path=os.path.join(current_path,"visualize/xplanes")

time_start=0
time_end=0
time_step=1
index_x=1
interpol=0 # 1 for interpol planes

var=["r","u","v","w","e"] #,"vortx","vorty","vortz"

x_scale=1
y_scale=1
z_scale=1

## functions

def getfluc(name,jmax,kmax,timestamp):
    c = 0
    data_time = np.zeros((jmax,kmax),precision)
    data_total = np.zeros((len(timestamp),jmax,kmax),precision)
    data_mean = np.zeros((jmax,kmax),precision)
    time_index=0
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_time=np.reshape(data, (jmax, kmax))
        data_total[time_index,:,:]=data_time
        time_index=time_index+1
    data_mean=np.mean(data_total,axis=0)
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_time=np.reshape(data, (jmax, kmax))
        data_fluc=data_time-data_mean
        data_fluc_reshape=np.reshape(data_fluc, (jmax*kmax))
        data_fluc_reshape.tofile("{0}.fluc.{1:07d}.bin".format(name,t)) 

precision = 'double'

if interpol==1:
    print('Loading interpolated planes')
    char_x='xpl_aI'
    x = np.fromfile('planes/x_I.bin', dtype=precision)
    y = np.fromfile('planes/y_I.bin', dtype=precision)
    z = np.fromfile('planes/z_I.bin', dtype=precision)
else:
    print('Loading current planes')
    char_x='xpl'
    x = np.fromfile('planes/x.bin', dtype=precision)
    y = np.fromfile('planes/y.bin', dtype=precision)
    z = np.fromfile('planes/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)

print(imax,jmax,kmax)

timestamps = np.arange(time_start,time_end+1,time_step)
print(timestamps)
for j in range(0,len(var)):
    getfluc(('planes/' + str(char_x) + '.' + str(index_x) + '.' + str(var[j])), jmax, kmax, timestamps)

datanames_var= ["" for j in range(len(var))]
datanames_fluc= ["" for j in range(len(var))]

for j in range(0,len(var)):
    datanames_var[j] = (str(char_x) + '.' + str(index_x) + "." + str(var[j]))
    datanames_fluc[j] = (str(char_x) + '.' + str(index_x) + "." + str(var[j]) + '.fluc')

ly = y[-1]+y[1]-y[0]

if interpol==1:
    for i in range(1, 5):
        writexmf("visualize/xplanes/xplanes_aI.{}.xmf".format(i), precision,  \
                [x[25]]*x_scale, (y+(i-1)*ly)*y_scale, z*z_scale,\
                timestamp = timestamps, dt = 1.0,\
                dataNames = datanames_var)
else:
    for i in range(1, 5):
        writexmf("visualize/xplanes/xplanes.{}.xmf".format(i), precision,  \
                [x[25]]*x_scale, (y+(i-1)*ly)*y_scale, z*z_scale,\
                timestamp = timestamps, dt = 1.0,\
                dataNames = datanames_var)

if interpol==1:
    writexmf("visualize/xplanes/xplanes_aI.fluc.xmf", precision,  \
            [x[25]]*x_scale, y*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_fluc)
else:
    writexmf("visualize/xplanes/xplanes.fluc.xmf", precision,  \
            [x[25]]*x_scale, y*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_fluc)



for i in range(0, len (timestamps)):
    print('timestamp={0:07d}'.format(timestamps[i]))
    timestamps_print='{0:07d}'.format(timestamps[i])
    for j in range(0,len(var)):
        src_path=os.path.join(current_path,"planes/" + str(char_x) + "."  + str(index_x) + "." + str(var[j]) + "." + str(timestamps_print) + ".bin")
        shutil.copy(src_path, dst_path) 
        src_path=os.path.join(current_path,"planes/" + str(char_x) + "."  + str(index_x) + "." + str(var[j]) +  ".fluc." + str(timestamps_print) + ".bin")
        shutil.copy(src_path, dst_path)

print('Planes copied')
