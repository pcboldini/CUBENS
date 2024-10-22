import shutil
import os
import numpy as np
from writexmf import writexmf

current_path=os.getcwd()
dst_path=os.path.join(current_path,"xplanes")

time_start=5
time_end=10
time_step=1
index_x=45
interpol=0 # 1 for interpol planes

var=["r","u","v","w","e"] #,"vortx","vorty","vortz"

x_scale=1
y_scale=1
z_scale=1

precision = 'double'

if interpol==1:
    print('Loading interpolated planes')
    char_x='xpl_aI'
    x = np.fromfile('../../output/planes/x_I.bin', dtype=precision)
    y = np.fromfile('../../output/planes/y_I.bin', dtype=precision)
    z = np.fromfile('../../output/planes/z_I.bin', dtype=precision)
else:
    print('Loading current planes')
    char_x='xpl'
    x = np.fromfile('../../output/planes/x.bin', dtype=precision)
    y = np.fromfile('../../output/planes/y.bin', dtype=precision)
    z = np.fromfile('../../output/planes/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)
print(imax,jmax,kmax)

slice_flag='off'
knew1=500
knew2=1000
jnew1=0
jnew2=jmax

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

## code

jmaxnew=jnew2-jnew1
kmaxnew=knew2-knew1

timestamps = np.arange(time_start,time_end+1,time_step)
print(timestamps)

for i in range(0, len (timestamps)):
    print('timestamp={0:07d}'.format(timestamps[i]))
    timestamps_print='{0:07d}'.format(timestamps[i])

if slice_flag=='on':
    for j in range(0,len(var)):
        for i in range(0, len (timestamps)):
            data=np.fromfile(os.path.join(current_path,"../../output/planes/" + str(char_x) + "." + str(index_x) + "." + str(var[j]) + "." + '{0:07d}'.format(timestamps[i]) + ".bin"))
            data_reshape=np.reshape(data, (kmax, jmax, 1))
            data_slice= data_reshape[knew1:knew2, jnew1:jnew2, :]
            data_slice_reshape = np.reshape(data_slice, (kmaxnew*1*jmaxnew))
            data_slice_reshape.tofile(os.path.join(current_path,"../../output/planes/" + str(char_x) + "." + str(index_x) + "." + str(var[j]) + ".slice." + '{0:07d}'.format(timestamps[i]) + ".bin"))

if slice_flag=='on':
    for j in range(0,len(var)):
        getfluc(('../../output/planes/' + str(char_x) + '.' + str(index_x) + '.'+ str(var[j]) + '.slice'), jmaxnew, kmaxnew, timestamps)
else:
    for j in range(0,len(var)):
        getfluc(('../../output/planes/' + str(char_x) + '.' + str(index_x) + '.' + str(var[j])), jmax, kmax, timestamps)


datanames_var= ["" for j in range(len(var))]
datanames_fluc= ["" for j in range(len(var))]

for j in range(0,len(var)):
    if slice_flag=='on':
        datanames_var[j]  =(str(char_x) + '.' + str(index_x) + '.' + str(var[j]) + '.slice')
    else:
        datanames_var[j]  =(str(char_x) + '.' + str(index_x) + '.' + str(var[j]))
    if slice_flag=='on':
        datanames_fluc[j] =(str(char_x) + '.' + str(index_x) + '.' + str(var[j]) + '.slice.fluc')
    else:
        datanames_fluc[j] =(str(char_x) + '.' + str(index_x) + '.' + str(var[j]) + '.fluc')

if slice_flag=='on':
    y=y[jnew1:jnew2]
    z=z[knew1:knew2]

ly = y[-1]+y[1]-y[0]

if interpol==1:
    for i in range(1, 5):
        writexmf("xplanes/xplanes_aI.{}.xmf".format(i), precision,  \
                [x[25]]*x_scale, (y+(i-1)*ly)*y_scale, z*z_scale,\
                timestamp = timestamps, dt = 1.0,\
                dataNames = datanames_var)
else:
    for i in range(1, 5):
        writexmf("xplanes/xplanes.{}.xmf".format(i), precision,  \
                [x[25]]*x_scale, (y+(i-1)*ly)*y_scale, z*z_scale,\
                timestamp = timestamps, dt = 1.0,\
                dataNames = datanames_var)

if interpol==1:
    writexmf("xplanes/xplanes_aI.fluc.xmf", precision,  \
            [x[25]]*x_scale, y*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_fluc)
else:
    writexmf("xplanes/xplanes.fluc.xmf", precision,  \
            [x[25]]*x_scale, y*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_fluc)


for j in range(0,len(var)):
    for i in range(0, len (timestamps)):
        if slice_flag=='on':
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_x) + "." + str(index_x) + "." + str(var[j]) + ".slice." + '{0:07d}'.format(timestamps[i]) + ".bin")
            shutil.copy(src_path, dst_path)
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_x) + "." + str(index_x) + "." + str(var[j]) +  ".slice.fluc." + '{0:07d}'.format(timestamps[i]) + ".bin")
            shutil.copy(src_path, dst_path)  
        else:
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_x) + "." + str(index_x) + "." + str(var[j]) + "." + '{0:07d}'.format(timestamps[i]) + ".bin")
            shutil.copy(src_path, dst_path)
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_x) + "." + str(index_x) + "." + str(var[j]) + ".fluc." + '{0:07d}'.format(timestamps[i]) + ".bin")
            shutil.copy(src_path, dst_path)     

print('Planes copied')
