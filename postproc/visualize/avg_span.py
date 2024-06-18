import shutil
import os
import numpy as np
from writexmf_one import writexmf_one

current_path=os.getcwd()
dst_path=os.path.join(current_path,"visualize/fields")

var=["r","T","u","v","w","q1","T13"]

x_scale=1
y_scale=1
z_scale=0.1

precision = 'double'
x = np.fromfile('../output/planes/x.bin', dtype=precision)
y = np.fromfile('../output/planes/y.bin', dtype=precision)
z = np.fromfile('../output/planes/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)

print(imax,jmax,kmax)

datanames_var= ["" for j in range(len(var))]

for j in range(0,len(var)):
    datanames_var[j] =('Yave_' + str(var[j]))

writexmf_one("visualize/fields/avg_span.xmf", precision,  \
        x*x_scale, [y[0]]*y_scale, z*z_scale,\
        timestamp = [0], dt = 0.0, \
        dataNames = datanames_var)

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


for j in range(0,len(var)):
    src_path=os.path.join(current_path,"results/Yave_" + str(var[j])  + ".bin")
    shutil.copy(src_path, dst_path)

print('Planes copied')

