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

path = '//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc'

z = np.array(netCDF4.Dataset(path).variables.get('lev'))
x = np.array(netCDF4.Dataset(path).variables.get('lon'))
delta_x = np.diff(x)
delta_x = np.insert(delta_x, 0, x[0])
u = np.array(netCDF4.Dataset(path).variables.get('u'))[:, :, 0, :]
v = np.array(netCDF4.Dataset(path).variables.get('v'))[:, :, 0, :]

bottom = np.where(z > 10)[0][0]
original_zdim = z.shape[0]
z2d = np.tile(z, (x.shape[0], 1)).T
z3d = np.tile(z, (u.shape[0], np.shape(z2d)[1], 1))
z = z3d.swapaxes(1, 2)
z = z[:, bottom:, :]
x2d = np.tile(x, (original_zdim, 1))
x3d = np.tile(x, (u.shape[0], np.shape(x2d)[0], 1))
x = x3d
delta_x2d = np.tile(delta_x, (original_zdim, 1))
delta_x3d = np.tile(delta_x, (u.shape[0], np.shape(x2d)[0], 1))
delta_x = delta_x3d
delta_x = delta_x[:, bottom:, :]
x = x[:, bottom:, :]
u = u[:, bottom:, :]
v = v[:, bottom:, :]

R = np.array(np.divide(np.multiply(np.sqrt(np.add(np.power(u, 2), np.power(v, 2))),
                                   np.multiply(delta_x, 1000)), 1.5 * (10 ** -5)))
R[R == np.inf] = np.nan
R_avg = np.nanmean(R, axis=0)

plt.imshow(R_avg, aspect='auto', cmap='nipy_spectral', vmin=0,
           extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average Reynolds Number')
plt.colorbar(label=r'Reynolds Number (unitless)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/rey_avg.png')
plt.close('all')

for i in range(R.shape[0]):
    plt.imshow(R[i, :, :], aspect='auto', cmap='nipy_spectral', vmin=0,
               extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'Reynolds Number at {i} Timestep')
    plt.colorbar(label=r'Reynolds Number (unitless)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/rey_hourly/rey_{i}.png')
    plt.close('all')
