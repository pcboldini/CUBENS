###############################################################################
# HEAD

import numpy as np
from writexmf import writexmf
import matplotlib
import matplotlib.pyplot as plt


###############################################################################
# PARAMETERS
path='planes'
end_time=1000000
delta_time=10000


###############################################################################
# FLUCTUATIONS

def getfluc(name,npts,timestamp):
    avg=np.zeros(npts,precision)
    
    avg = np.fromfile("{0}.{1:07d}.bin".format(name,0), dtype=precision)
    for t in timestamp:
        data = avg - np.fromfile("{0}.{1:07d}.bin".format(name,t), dtype=precision) 
        data.tofile("{0}.fluc.{1:07d}.bin".format(name,t))

###############################################################################


###############################################################################
# XMF for paraview
        
precision = 'double'
x = np.fromfile('../../' + path + '/x.bin', dtype=precision)
y = np.fromfile('../../' + path + '/y.bin', dtype=precision)
z = np.fromfile('../../' + path + '/z.bin', dtype=precision)

imax = np.size(x)
jmax = np.size(y)
kmax = np.size(z)

print('Number of points')
print('z-direction: {}'.format(kmax))
print('x-direction: {}'.format(imax))
print('y-direction: {}'.format(jmax))

timestamps = np.arange(0,end_time,delta_time)
number_time=np.size(timestamps)
print('\nTime extractions')
print('End time: {}'.format(end_time))
print('Delta time: {}'.format(delta_time))
print('Count: {}'.format(number_time))

# Total variables
writexmf("../../results/yplanes.xmf", precision,  \
        x, [y[0]], z*0.1,\
        timestamp = timestamps, dt = 1.0,\
        dataNames = ['../' + path + '/ypl.r',\
                     '../' + path + '/ypl.p',\
                     '../' + path + '/ypl.t',\
                     '../' + path + '/ypl.e',\
                     '../' + path + '/ypl.u',\
                     '../' + path + '/ypl.w'])

# Fluctuations    
#getfluc('../../planes/ypl.r', kmax*jmax, timestamps)
#getfluc('../../planes/ypl.p', kmax*jmax, timestamps)
#getfluc('../../planes/ypl.t', kmax*jmax, timestamps)
#getfluc('../../planes/ypl.u', kmax*jmax, timestamps)
#getfluc('../../planes/ypl.w', kmax*jmax, timestamps)
    
#writexmf("../../results/yplanes.fluc.xmf", precision,  \
#        x, [y[0]], z,\
#        timestamp = timestamps, dt = 1.0,\
#        dataNames = ['planes/ypl.r.fluc',\
#                     'planes/ypl.p.fluc',\
#                     'planes/ypl.t.fluc',\
#                     'planes/ypl.e.fluc',\
#                     'planes/ypl.u.fluc',\
#                     'planes/ypl.w.fluc'])
    


