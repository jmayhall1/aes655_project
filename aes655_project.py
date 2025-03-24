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
data_change = np.array(netCDF4.Dataset(path).variables.get('dqke'))[:, bottom:, 0, :]
data_adv = np.array(netCDF4.Dataset(path).variables.get('qke_adv'))[:, bottom:, 0, :]
data_ke = np.array(netCDF4.Dataset(path).variables.get('qke'))[:, bottom:, 0, :]
data_vert = np.array(netCDF4.Dataset(path).variables.get('qwt'))[:, bottom:, 0, :]

data_avg_shear = np.nanmean(data_shear, axis=0)
data_avg_buoy = np.nanmean(data_buoy, axis=0)
data_avg_diss = np.nanmean(data_diss, axis=0)
data_avg_change = np.nanmean(data_change, axis=0)
data_avg_adv = np.nanmean(data_adv, axis=0)
data_avg_ke = np.nanmean(data_ke, axis=0)
data_avg_vert = np.nanmean(data_vert, axis=0)

plt.imshow(data_avg_shear, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE Shear Production')
plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/shear_avg.png')
plt.close('all')

plt.imshow(data_avg_buoy, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE Buoyancy Production')
plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/buoy_avg.png')
plt.close('all')

plt.imshow(data_avg_diss, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE Dissipation Production')
plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/diss_avg.png')
plt.close('all')

plt.imshow(data_avg_change, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE Change')
plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/change_avg.png')
plt.close('all')

plt.imshow(data_avg_adv, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE Advection')
plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/adv_avg.png')
plt.close('all')

plt.imshow(data_avg_vert, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE Vertical Transport')
plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/vert_avg.png')
plt.close('all')

plt.imshow(data_avg_ke, aspect='auto', cmap='rainbow',
           extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel(r'Range (km)')
plt.title(f'Average TKE')
plt.colorbar(label=r'TKE ($\frac{m^2}{s^2}$)')
# plt.show()
plt.savefig(f'//uahdata/rstor/aes655_project/ke_avg.png')
plt.close('all')

for i in range(data_shear.shape[0]):
    plt.imshow(data_shear[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE Shear Production at {i} Timestep')
    plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/shear_hourly/shear_{i}.png')
    plt.close('all')

    plt.imshow(data_buoy[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE Buoyancy Production at {i} Timestep')
    plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/buoy_hourly/buoy_{i}.png')
    plt.close('all')

    plt.imshow(data_diss[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE Dissipation Production at {i} Timestep')
    plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/diss_hourly/diss_{i}.png')
    plt.close('all')

    plt.imshow(data_change[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE Change at {i} Timestep')
    plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/change_hourly/change_{i}.png')
    plt.close('all')

    plt.imshow(data_adv[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE Advection at {i} Timestep')
    plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/adv_hourly/adv_{i}.png')
    plt.close('all')

    plt.imshow(data_vert[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE Vertical Transport at {i} Timestep')
    plt.colorbar(label=r'TKE per Second ($\frac{m^2}{s^3}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/vert_hourly/vert_{i}.png')
    plt.close('all')

    plt.imshow(data_ke[i, :, :], aspect='auto', cmap='rainbow',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel(r'Range (km)')
    plt.title(f'TKE at {i} Timestep')
    plt.colorbar(label=r'TKE ($\frac{m^2}{s^2}$)')
    # plt.show()
    plt.savefig(f'//uahdata/rstor/aes655_project/ke_hourly/ke_{i}.png')
    plt.close('all')