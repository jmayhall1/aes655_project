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

path = '//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc'

# data = netCDF4.Dataset(path).variables
z = np.array(netCDF4.Dataset(path).variables.get('lev'))
x = np.array(netCDF4.Dataset(path).variables.get('lon'))
bottom = np.where(z > 10)[0][0]
z = z[bottom:]
data_shear = np.array(netCDF4.Dataset(path).variables.get('qshear'))[:, bottom:, 0, :]
data_buoy = np.array(netCDF4.Dataset(path).variables.get('qbuoy'))[:, bottom:, 0, :]
data_diss = np.array(netCDF4.Dataset(path).variables.get('qdiss'))[:, bottom:, 0, :]

buoy_shear = np.abs(np.divide(data_buoy, data_shear)) * 100
buoy_shear[buoy_shear == np.inf] = np.nan
buoy_avg_shear = np.abs(np.nanmean(buoy_shear, axis=0))
diss_shear = np.abs(np.divide(data_diss, data_shear)) * 100
diss_shear[diss_shear == np.inf] = np.nan
diss_avg_shear = np.abs(np.nanmean(diss_shear, axis=0))
buoydiss_shear = np.abs(np.divide(np.add(data_buoy, data_diss), data_shear)) * 100
buoydiss_shear[buoydiss_shear == np.inf] = np.nan
buoydiss_avg_shear = np.abs(np.nanmean(buoydiss_shear, axis=0))

plt.imshow(buoy_avg_shear, aspect='auto', cmap='rainbow', vmin=0, vmax=110,
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average Buoyancy Production %')
plt.colorbar(label=r'%')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/buoy_avg_percent.png')
plt.close('all')

plt.imshow(diss_avg_shear, aspect='auto', cmap='rainbow', vmin=0, vmax=110,
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average Dissipation %')
plt.colorbar(label=r'%')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/diss_avg_percent.png')
plt.close('all')

plt.imshow(buoydiss_avg_shear, aspect='auto', cmap='rainbow', vmin=0, vmax=110,
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average Buoyancy Production and Dissipation %')
plt.colorbar(label=r'%')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/buoydiss_avg_percent.png')
plt.close('all')

for i in range(data_shear.shape[0]):
    plt.imshow(buoy_shear[i, :, :], aspect='auto', cmap='rainbow', vmin=0, vmax=110,
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'Buoyancy Production % at {i} Timestep')
    plt.colorbar(label=r'%')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/buoy_percent_hourly/buoy_percent_{i}.png')
    plt.close('all')

    plt.imshow(diss_shear[i, :, :], aspect='auto', cmap='rainbow', vmin=0, vmax=110,
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'Dissipation % at {i} Timestep')
    plt.colorbar(label=r'%')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/diss_percent_hourly/diss_percent_{i}.png')
    plt.close('all')

    plt.imshow(buoydiss_shear[i, :, :], aspect='auto', cmap='rainbow', vmin=0, vmax=110,
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'Buoyancy Production and Dissipation % at {i} Timestep')
    plt.colorbar(label=r'%')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/buoydiss_percent_hourly/buoydiss_percent_{i}.png')
    plt.close('all')
