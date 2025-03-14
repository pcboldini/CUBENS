import shutil
import os
import numpy as np
from writexmf import writexmf

# Remove all files inside the yplanes folder
dst_path = os.path.join(os.getcwd(), "yplanes")
for filename in os.listdir(dst_path):
    file_path = os.path.join(dst_path, filename)
    if os.path.isfile(file_path):
        os.remove(file_path)

current_path=os.getcwd()

time_start=5
time_end=10
time_step=1
index_y=1
interpol=0 # 1 for interpol planes

var=["r","w"]

x_scale=1
y_scale=1
z_scale=1

precision = 'double'

if interpol==1:
    print('Loading interpolated planes')
    char_y='ypl_aI'
    x = np.fromfile('../../output/planes/x_I.bin', dtype=precision)
    y = np.fromfile('../../output/planes/y_I.bin', dtype=precision)
    z = np.fromfile('../../output/planes/z_I.bin', dtype=precision)
else:
    print('Loading current planes')
    char_y='ypl'
    x = np.fromfile('../../output/planes/x.bin', dtype=precision)
    y = np.fromfile('../../output/planes/y.bin', dtype=precision)
    z = np.fromfile('../../output/planes/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)
print(imax,jmax,kmax)

slice_flag='on'
inew1=0
inew2=100
knew1=1
knew2=500

## functions

def getfluc(name,imax,kmax,timestamp):
    c = 0
    data_time = np.zeros((imax,kmax),precision)
    data_total = np.zeros((len(timestamp),imax,kmax),precision)
    data_mean = np.zeros((imax,kmax),precision)
    time_index=0
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_time=np.reshape(data, (imax, kmax))
        data_total[time_index,:,:]=data_time
        time_index=time_index+1
    data_mean=np.mean(data_total,axis=0)
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_time=np.reshape(data, (imax, kmax))
        data_fluc=data_time-data_mean
        data_fluc_reshape=np.reshape(data_fluc, (imax*kmax))
        data_fluc_reshape.tofile("{0}.fluc.{1:07d}.bin".format(name,t)) 


## code

imaxnew=inew2-inew1
kmaxnew=knew2-knew1

timestamps = np.arange(time_start,time_end+1,time_step)
print(timestamps)

for i in range(0, len (timestamps)):
    print('timestamp={0:07d}'.format(timestamps[i]))
    timestamps_print='{0:07d}'.format(timestamps[i])

if slice_flag=='on':
    for j in range(0,len(var)):
        for i in range(0, len (timestamps)):
            data=np.fromfile(os.path.join(current_path,"../../output/planes/" + str(char_y) + "." + str(index_y) + "." + str(var[j]) + "." + '{0:07d}'.format(timestamps[i]) + ".bin"))
            data_reshape=np.reshape(data, (kmax, 1, imax))
            data_slice= data_reshape[knew1:knew2, :, inew1:inew2]
            data_slice_reshape = np.reshape(data_slice, (kmaxnew*1*imaxnew))
            data_slice_reshape.tofile(os.path.join(current_path,"../../output/planes/" + str(char_y) + "." + str(index_y) + "." + str(var[j]) + ".slice." + '{0:07d}'.format(timestamps[i]) + ".bin"))

if slice_flag=='on':
    for j in range(0,len(var)):
        getfluc(('../../output/planes/' + str(char_y) + '.' + str(index_y) + '.'+ str(var[j]) + '.slice'), imaxnew, kmaxnew, timestamps)
else:
    for j in range(0,len(var)):
        getfluc(('../../output/planes/' + str(char_y) + '.' + str(index_y) + '.'+ str(var[j])), imax, kmax, timestamps)

datanames_var= ["" for j in range(len(var))]
datanames_fluc= ["" for j in range(len(var))]

for j in range(0,len(var)):
    if slice_flag=='on':
        datanames_var[j]  =(str(char_y) + '.' + str(index_y) + '.' + str(var[j]) + '.slice')
    else:
        datanames_var[j]  =(str(char_y) + '.' + str(index_y) + '.' + str(var[j]))
    if slice_flag=='on':
        datanames_fluc[j] =(str(char_y) + '.' + str(index_y) + '.' + str(var[j]) + '.slice.fluc')
    else:
        datanames_fluc[j] =(str(char_y) + '.' + str(index_y) + '.' + str(var[j]) + '.fluc')

xa = 20

if slice_flag=='on':
    x=x[inew1:inew2]
    z=z[knew1:knew2]

if interpol==1:
    writexmf("yplanes/yplanes_aI.xmf", precision,  \
            x*x_scale-xa, [y[0]]*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_var)
else:
    writexmf("yplanes/yplanes.xmf", precision,  \
            x*x_scale, [y[0]]*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_var)

if interpol==1:
    writexmf("yplanes/yplanes_aI.fluc.xmf", precision,  \
            x*x_scale, [y[0]]*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_fluc)
else:
    writexmf("yplanes/yplanes.fluc.xmf", precision,  \
            x*x_scale, [y[0]]*y_scale, z*z_scale,\
            timestamp = timestamps, dt = 1.0,\
            dataNames = datanames_fluc)


for j in range(0,len(var)):
    for i in range(0, len (timestamps)):
        if slice_flag=='on':
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_y) + "." + str(index_y) + "." + str(var[j]) + ".slice." + '{0:07d}'.format(timestamps[i])+ ".bin")
            shutil.move(src_path, dst_path)
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_y) + "." + str(index_y) + "." + str(var[j]) +  ".slice.fluc." + '{0:07d}'.format(timestamps[i]) + ".bin")
            shutil.move(src_path, dst_path)  
        else:
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_y) + "." + str(index_y) + "." + str(var[j]) +  "." + '{0:07d}'.format(timestamps[i])+ ".bin")
            shutil.move(src_path, dst_path)
            src_path=os.path.join(current_path,"../../output/planes/" + str(char_y) + "." + str(index_y) + "." + str(var[j]) +  ".fluc." + '{0:07d}'.format(timestamps[i])+ ".bin")
            shutil.move(src_path, dst_path)       


print('Planes copied')
