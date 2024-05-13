import shutil
import os
import numpy as np
from writexmf import writexmf

current_path=os.getcwd()
dst_path=os.path.join(current_path,"visualize/fields")

time_start=0
time_end=100
time_step=10

var1=["w"] # r,u,v,w,e
var2="qvort"

x_scale=1
y_scale=1
z_scale=1

precision = 'double'

x = np.fromfile('planes/x.bin', dtype=precision)
y = np.fromfile('planes/y.bin', dtype=precision)
z = np.fromfile('planes/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)

print(imax,jmax,kmax)

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
            reshape=np.reshape(data[index_var[tt]*ntot:(index_var[tt]+1)*ntot], (imax, jmax, kmax))
            reshape.tofile(os.path.join("{0}."+ str(var[tt])+ ".{1:07d}.bin").format(name,t))
   
timestamps = np.arange(time_start,time_end,time_step)
array_var=["r","u","v","w","e"]
index_var1 = np.zeros(len(var1),'int')
for i in range(0,len(array_var)):
    for j in range(0,len(var1)):
    	if array_var[i] == var1[j]:
           index_var1[j]=i

getfields(('../restart/' + 'ruvwe' ), var1, index_var1, imax, jmax, kmax, timestamps)

datanames_var= ["" for j in range(len(var1))]

for j in range(0,len(var1)):
    datanames_var[j] =('ruvwe' + '.' + str(var1[j]))

datanames_var.append(var2)

print(datanames_var)

ly = y[-1]+y[1]-y[0]

for i in range(1, 5):
    writexmf("visualize/fields/3dvars.{}.xmf".format(i), precision,  \
                x*x_scale, (y+(i-1)*ly)*y_scale, z*z_scale,\
                timestamp = timestamps, dt = 1.0,\
                dataNames = datanames_var)


for i in range(0, len (timestamps)):
    print('timestamp={0:07d}'.format(timestamps[i]))
    timestamps_print='{0:07d}'.format(timestamps[i])
    for j in range(0,len(var1)):
        src_path=os.path.join(current_path,"../restart/ruvwe." + str(var1[j]) + "." + str(timestamps_print) + ".bin")
        shutil.copy(src_path, dst_path)
    src_path=os.path.join(current_path,"results/qvort." + str(timestamps_print) + ".bin")
    shutil.copy(src_path, dst_path)

print('Planes copied')
