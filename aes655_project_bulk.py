# coding=utf-8
"""
@author: John Mark Mayhall
"""
import glob
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
from mayavi import mlab
import metpy.calc as mpcalc

# theta_prime = mpcalc.first_derivative(vpt, x=z)
# u_prime = mpcalc.first_derivative(uz, x=z)
# v_prime = mpcalc.first_derivative(vz, x=z)

path = '//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc'

# data = netCDF4.Dataset(path).variables
z = np.array(netCDF4.Dataset(path).variables.get('lev'))
delta_z = np.diff(z)
delta_z = np.insert(delta_z, 0, z[0])
x = np.array(netCDF4.Dataset(path).variables.get('lon'))
u = np.array(netCDF4.Dataset(path).variables.get('u'))[:, :, 0, :]
v = np.array(netCDF4.Dataset(path).variables.get('v'))[:, :, 0, :]
thv = np.array(netCDF4.Dataset(path).variables.get('thv'))[:, :, 0, :]
delta_u = mpcalc.first_derivative(u, axis=1, x=z)
delta_v = mpcalc.first_derivative(v, axis=1, x=z)
delta_thv = mpcalc.first_derivative(thv, axis=1, x=z)

bottom = np.where(z > 10)[0][0]
z2d = np.tile(z, (x.shape[0], 1)).T
z3d = np.tile(z, (u.shape[0], np.shape(z2d)[1], 1))
z = z3d.swapaxes(1, 2)
z = np.multiply(z[:, bottom:, :], 1000)
delta_z2d = np.tile(delta_z, (x.shape[0], 1)).T
delta_z3d = np.tile(delta_z, (u.shape[0], np.shape(z2d)[1], 1))
delta_z = delta_z3d.swapaxes(1, 2)
delta_z = np.multiply(delta_z[:, bottom:, :], 1000)
u = u[:, bottom:, :]
v = v[:, bottom:, :]
thv = thv[:, bottom:, :]
delta_u = delta_u[:, bottom:, :]
delta_v = delta_v[:, bottom:, :]
delta_thv = delta_thv[:, bottom:, :]

R = np.array(9.81 * delta_z * delta_thv / (thv * (delta_u ** 2 + delta_v ** 2)))
R[R == np.inf] = np.nan
R_avg = np.nanmean(R, axis=0)

plt.imshow(R_avg, aspect='auto', cmap='rainbow', vmin=-1, vmax=10,
           extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Longitude ($^\circ$)')
plt.title(f'Average Richardson Number')
plt.colorbar(label=r'$R_b$ (unitless)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/Rb_avg.png')
plt.close('all')

time = []
minimum_value = []
minimum_loc = np.zeros((30, 100))
loc025 = np.zeros((30, 100))
loc1 = np.zeros((30, 100))
for i in range(R.shape[0]):
    time.append(i)
    minimum_value.append(np.nanmin(R[i, :, :]))
    current_array = R[i, :, :] == np.nanmin(R[i, :, :])
    array025 = R[i, :, :] < 0.25
    array1 = R[i, :, :] < 1
    minimum_loc = np.add(minimum_loc, current_array.astype(int))
    loc025 = np.add(loc025, array025.astype(int))
    loc1 = np.add(loc1, array1.astype(int))
    plt.imshow(R[i, :, :], aspect='auto', cmap='rainbow', vmin=-1, vmax=1,
               extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Longitude ($^\circ$)')
    plt.title(f'$R_b$ at {i} Timestep')
    plt.colorbar(label=r'$R_b$ (unitless)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/Rb_hourly/Rb_{i}.png')
    plt.close('all')

plt.plot(time, minimum_value)
plt.ylabel('$R_b$')
plt.xlabel(r'Time Step')
plt.ylim(0, 2)
plt.title(f' Minimum $R_b$ vs Time')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/Rb_min.png')
plt.close('all')

plt.imshow(minimum_loc, aspect='auto', cmap='rainbow', vmin=0,
           extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Longitude ($^\circ$)')
plt.title(f'Location of Minimum $R_b$')
plt.colorbar(label=r'# of Minimum $R_b$ Occurrences')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/minRb_loc.png')
plt.close('all')

plt.imshow(loc025, aspect='auto', cmap='rainbow', vmin=0,
           extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Longitude ($^\circ$)')
plt.title(f'Location of $R_b$ < 0.25')
plt.colorbar(label=r'# of $R_b$ < 0.25 Occurrences')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/025Rb_loc.png')
plt.close('all')

plt.imshow(loc1, aspect='auto', cmap='rainbow', vmin=0,
           extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Longitude ($^\circ$)')
plt.title(f'Location of $R_b$ < 1')
plt.colorbar(label=r'# of $R_b$ < 1 Occurrences')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/1Rb_loc.png')
plt.close('all')
