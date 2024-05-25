import shutil
import os
import numpy as np
from writexmf import writexmf

current_path=os.getcwd()
dst_path=os.path.join(current_path,"visualize/zplanes")

time_start=0
time_end=0
time_step=1
index_z=500

var=["r","p","t","u","v","w"]
#var = ["w"]

x_scale=1
y_scale=1
z_scale=1

## functions

def getfluc(name,jmax,imax,timestamp):
    c = 0
    data_time = np.zeros((jmax,imax),precision)
    data_total = np.zeros((len(timestamp),jmax,imax),precision)
    data_mean = np.zeros((jmax,imax),precision)
    time_index=0
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_time=np.reshape(data, (jmax,imax))
        data_total[time_index,:,:]=data_time
        time_index=time_index+1
    data_mean=np.mean(data_total,axis=0)
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_time=np.reshape(data, (jmax,imax))
        data_fluc=data_time-data_mean
        data_fluc_reshape=np.reshape(data_fluc, (imax*jmax))
        data_fluc_reshape.tofile("{0}.fluc.{1:07d}.bin".format(name,t)) 

precision = 'double'
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
    getfluc(('planes/zpl.' + str(index_z) + '.'+ str(var[j])), jmax, imax, timestamps)

datanames_var= ["" for j in range(len(var))]
datanames_fluc= ["" for j in range(len(var))]

for j in range(0,len(var)):
    datanames_var[j] =('zpl.'  + str(index_z) + '.' + str(var[j]))
    datanames_fluc[j] =('zpl.' + str(index_z) + '.' + str(var[j]) + '.fluc')

writexmf("visualize/zplanes/zplanes.xmf", precision,  \
        x*x_scale, y*y_scale, [z[0]]*z_scale,\
        timestamp = timestamps, dt = 1.0,\
        dataNames = datanames_var)

writexmf("visualize/zplanes/zplanes.fluc.xmf", precision,  \
        x*x_scale, y*y_scale, [z[0]]*z_scale,\
        timestamp = timestamps, dt = 1.0,\
        dataNames = datanames_fluc)

#xa = 22
#for i in range(0, 5):
#    writexmf("visualize_fields/yplane2.{}.xmf".format(i), precision, \
#             x-i*xa, [y[0]], z*1.0, \
#             timestamp = timestamps, dt = 1.0, \
#             dataNames = ['ypl.u.fluc',\
#                          'ypl.v.fluc',\
#                          'ypl.w.fluc',\
#                          'ypl.r.fluc',\
#                          'ypl.p.fluc'])

for i in range(0, len (timestamps)):
    print('timestamp={0:07d}'.format(timestamps[i]))
    timestamps_print='{0:07d}'.format(timestamps[i])

    for j in range(0,len(var)):
        src_path=os.path.join(current_path,"planes/zpl." + str(index_z) + "." + str(var[j]) + "." + str(timestamps_print) + ".bin")
        shutil.copy(src_path, dst_path)
        src_path=os.path.join(current_path,"planes/zpl." + str(index_z) + "." + str(var[j]) +  ".fluc." + str(timestamps_print) + ".bin")
        shutil.copy(src_path, dst_path)

print('Planes copied')
