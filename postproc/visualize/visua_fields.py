import shutil
import os
import numpy as np
from writexmf import writexmf

current_path=os.getcwd()
dst_path=os.path.join(current_path,"fields")

time_start=5
time_end=10
time_step=1

var1=["w"] # r,u,v,w,e
var2="qvort"

x_scale=1
y_scale=1
z_scale=1

precision = 'double'

x = np.fromfile('../../output/planes/x.bin', dtype=precision)
y = np.fromfile('../../output/planes/y.bin', dtype=precision)
z = np.fromfile('../../output/planes/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)
print(imax,jmax,kmax)

slice_flag='on'
inew1=0
inew2=200
knew1=0
knew2=2000

# Functions

def getfields(name,var,index_var,imax,jmax,kmax,timestamp):
    ntot=imax*jmax*kmax
    for t in timestamp:
        data=np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision)
        data_head=data[0:4]
        data=data[5:]
        data_r=data[0:ntot]
        data_u=data[ntot:2*ntot]
        data_v=data[2*ntot:3*ntot]
        data_w=data[3*ntot:4*ntot]
        data_e=data[4*ntot:5*ntot]
        for tt in range(0,len(index_var)):
            reshape=np.reshape(data[index_var[tt]*ntot:(index_var[tt]+1)*ntot], (kmax, jmax, imax))
            if slice_flag=='on':
                data_slice= reshape[knew1:knew2, :, inew1:inew2]
                data_slice_reshape = np.reshape(data_slice, (kmaxnew*jmax*imaxnew))
                data_slice_reshape.tofile(os.path.join("{0}."+ str(var[tt])+ ".slice.{1:07d}.bin").format(name,t))
            reshape.tofile(os.path.join("{0}."+ str(var[tt])+ ".{1:07d}.bin").format(name,t))

# Code

imaxnew=inew2-inew1
kmaxnew=knew2-knew1

if time_start==time_end:
    timestamps = np.arange(time_start,time_end+1,time_step)
else:
    timestamps = np.arange(time_start,time_end,time_step) 
    
for i in range(0, len (timestamps)):
    print('timestamp={0:07d}'.format(timestamps[i]))
    timestamps_print='{0:07d}'.format(timestamps[i])

array_var=["r","u","v","w","e"]
index_var1 = np.zeros(len(var1),'int')
for i in range(0,len(array_var)):
    for j in range(0,len(var1)):
    	if array_var[i] == var1[j]:
           index_var1[j]=i

getfields(('../../output/restart/' + 'ruvwe' ), var1, index_var1, imax, jmax, kmax, timestamps)

if slice_flag=='on':
    for i in range(0, len (timestamps)):
        data=np.fromfile(os.path.join(current_path,"../results/vort/qvort." + '{0:07d}'.format(timestamps[i]) + ".bin"))
        data_reshape=np.reshape(data, (kmax, jmax, imax))
        data_slice= data_reshape[knew1:knew2, :, inew1:inew2]
        data_slice_reshape = np.reshape(data_slice, (kmaxnew*jmax*imaxnew))
        data_slice_reshape.tofile(os.path.join(current_path,"../results/vort/qvort.slice." + '{0:07d}'.format(timestamps[i]) + ".bin"))

datanames_var= ["" for j in range(len(var1))]

for j in range(0,len(var1)):
    if slice_flag=='on':
        datanames_var[j] =('ruvwe' + '.' + str(var1[j]) + '.slice')
    else:
        datanames_var[j] =('ruvwe' + '.' + str(var1[j]))

if slice_flag=='on':
    datanames_var.append(var2+'.slice')
else:
    datanames_var.append(var2)

print(datanames_var)

if slice_flag=='on':
    x=x[inew1:inew2]
    z=z[knew1:knew2]

ly = y[-1]+y[1]-y[0]

for i in range(1, 5):
    writexmf("fields/3dvars.{}.xmf".format(i), precision,  \
                x*x_scale, (y+(i-1)*ly)*y_scale, z*z_scale,\
                timestamp = timestamps, dt = 1.0,\
                dataNames = datanames_var)

for i in range(0, len (timestamps)):
    for j in range(0,len(var1)):
        if slice_flag=='on':
            src_path=os.path.join(current_path,"../../output/restart/ruvwe." + str(var1[j]) + ".slice." + '{0:07d}'.format(timestamps[i])+ ".bin")
            shutil.copy(src_path, dst_path)
        else:
            src_path=os.path.join(current_path,"../../output/restart/ruvwe." + str(var1[j]) +  "." + '{0:07d}'.format(timestamps[i])+ ".bin")
            shutil.copy(src_path, dst_path)
    if slice_flag=='on':
        src_path=os.path.join(current_path,"../results/vort/qvort.slice."+ '{0:07d}'.format(timestamps[i])+ ".bin")
    else:
        src_path=os.path.join(current_path,"../results/vort/qvort." + '{0:07d}'.format(timestamps[i])+ ".bin")   
    shutil.copy(src_path, dst_path)

print('Planes copied')
