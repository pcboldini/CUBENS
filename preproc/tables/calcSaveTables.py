# -
#
# SPDX-FileCopyrightText: Copyright (c) 2024 Pietro Carlo Boldini, Rene Pecnik and the CUBENS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -

from makePropertyTable import *



params = {'fluid' : 'REFPROP::CO2', 
          'Ec'    : 0.01, 
          'P_inf' : 85.0e5, 
          'T_inf' : 294.9593531839261, 
          'T_wall': 326.0077061506552,
          'Re'    : 2162.6838776652576}


inf = refValues(params)



RE_tab = RE_table(params['fluid'], 
                   Tmin = 280.0, Tmax = 340.0, 
                   Pmin = 78e5,  Pmax = 95e5, 
                   imax = 400, jmaxDiag=140, jOffset=60)

RE_tab.saveTable(inf)

fig, ax = plt.subplots(1,1, figsize=(14,10))
RE_tab.plotTable_RE(fig, ax, params['T_inf'], params['T_wall'], params['P_inf'], RE_tab.tem, 0)




# # Make $P-T$ table

imax = 55
jmax = 600

pre_1d = np.linspace(80e5 , 90e5,  imax)
tem_1d = np.linspace(270.0, 350.0, jmax)
pre,tem = np.meshgrid(pre_1d, tem_1d, indexing='ij')
rho = np.zeros((imax, jmax))

for i,p in enumerate(pre_1d):
    rho[i,:] = CP.PropsSI('D', "P", p, "T", tem_1d, params['fluid'])

( pre_1d / (inf.rho*inf.vel**2.0) ).tofile("PT_tab_1d_pre.bin")
( tem_1d / inf.tem                ).tofile("PT_tab_1d_tem.bin")
( rho    / inf.rho                ).tofile("PT_tab_2d_rho.bin")


fig, ax = plt.subplots(1,1,figsize=(14,6))
plt.contourf(tem, pre, rho, 20, cmap = "jet")
plt.colorbar()
plt.plot(tem, pre,     'k', lw=0.5)
plt.plot(tem.T, pre.T, 'k', lw=0.5)
plt.plot([params['T_inf'], params['T_wall']], params['P_inf']*np.ones(2), 'k-',lw=3, label="BL states")




# # Make $\rho-T$ table

imax = 300
jmax = 400

rmin = np.min(rho); print(rmin)
rmax = np.max(rho); print(rmax)

Tmin = 270.0
Tmax = 340.0

rho_1d = np.linspace(rmin, rmax, imax)
tem_1d = np.linspace(Tmin, Tmax, jmax)
rho,tem = np.meshgrid(rho_1d, tem_1d, indexing='ij')
ien = np.zeros((imax, jmax))
for i,r in enumerate(rho_1d):
    ien[i,:] = CP.PropsSI('U', "D", r, "T", tem_1d, params['fluid'])

( rho_1d / inf.rho    ).tofile("RT_tab_1d_rho.bin")
( tem_1d / inf.tem    ).tofile("RT_tab_1d_tem.bin")
( ien    / inf.vel**2 ).tofile("RT_tab_2d_ien.bin")


fig, ax = plt.subplots(1,1,figsize=(14,10))
plt.contourf(rho, tem, ien, 20, cmap="jet")
plt.plot(rho,   (tem),   'k', lw=0.5)
plt.plot(rho.T, (tem.T), 'k', lw=0.5)
plt.colorbar()

# reference conditions of boundary layer
T_bl = np.linspace( params['T_inf'], params['T_wall'], 100)
plt.plot(CP.PropsSI('D', 'P', params['P_inf'], 'T', T_bl, params['fluid']), T_bl,   
            'k-', lw=3, label="Boundary layer distribution")


plt.show()
