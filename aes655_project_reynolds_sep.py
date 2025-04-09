# coding=utf-8
"""
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np

if __name__ == "__main__":
    # Define file paths
    base_path = '//uahdata/rstor/aes655_project/'

    # Load NetCDF data efficiently
    ncfile = netCDF4.Dataset(os.path.join(base_path, 'cm1out_azimavg_s.nc'))
    base_path = '//uahdata/rstor/aes655_project/sep_by_intensity_phase/'
    z = np.asarray(ncfile.variables['lev'])
    x = np.asarray(ncfile.variables['lon'])
    mtime = (np.asarray(ncfile.variables['mtime']) / 3600)[:, 0, 0]
    pre_ri_time = np.where(mtime == 10)[0][0]
    post_ri_time = np.where(mtime == 71)[0][0]
    preri_u = np.asarray(ncfile.variables['u'])[:pre_ri_time, :, 0, :]
    preri_v = np.asarray(ncfile.variables['v'])[:pre_ri_time, :, 0, :]
    ri_u = np.asarray(ncfile.variables['u'])[pre_ri_time:post_ri_time, :, 0, :]
    ri_v = np.asarray(ncfile.variables['v'])[pre_ri_time:post_ri_time, :, 0, :]
    postri_u = np.asarray(ncfile.variables['u'])[post_ri_time:, :, 0, :]
    postri_v = np.asarray(ncfile.variables['v'])[post_ri_time:, :, 0, :]
    ncfile.close()  # Close file after loading data

    # Compute spatial step differences
    delta_x = np.diff(x, prepend=x[0])

    bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
    top = np.searchsorted(z, 20) + 1

    # Reshape arrays efficiently
    original_zdim, xdim = z.shape[0], x.shape[0]
    z3d = np.tile(z[:, np.newaxis], (1, xdim))[np.newaxis, :, :]
    x3d = np.tile(x[np.newaxis, :], (original_zdim, 1))
    delta_x3d = np.tile(delta_x[np.newaxis, :], (original_zdim, 1))

    preri_z3d = np.repeat(z3d, preri_u.shape[0], axis=0)[:, bottom: top, :]
    preri_x3d = np.repeat(x3d[np.newaxis, :, :], preri_u.shape[0], axis=0)[:, bottom: top, :]
    preri_delta_x3d = np.repeat(delta_x3d[np.newaxis, :, :], preri_u.shape[0], axis=0)[:, bottom: top, :]
    ri_z3d = np.repeat(z3d, ri_u.shape[0], axis=0)[:, bottom: top, :]
    ri_x3d = np.repeat(x3d[np.newaxis, :, :], ri_u.shape[0], axis=0)[:, bottom: top, :]
    ri_delta_x3d = np.repeat(delta_x3d[np.newaxis, :, :], ri_u.shape[0], axis=0)[:, bottom: top, :]
    postri_z3d = np.repeat(z3d, postri_u.shape[0], axis=0)[:, bottom: top, :]
    postri_x3d = np.repeat(x3d[np.newaxis, :, :], postri_u.shape[0], axis=0)[:, bottom: top, :]
    postri_delta_x3d = np.repeat(delta_x3d[np.newaxis, :, :], postri_u.shape[0], axis=0)[:, bottom: top, :]

    # Slice u and v
    preri_u, preri_v = preri_u[:, bottom: top, :], preri_v[:, bottom: top, :]
    ri_u, ri_v = ri_u[:, bottom: top, :], ri_v[:, bottom: top, :]
    postri_u, postri_v = postri_u[:, bottom: top, :], postri_v[:, bottom: top, :]

    # Compute Reynolds number
    speed_sq = preri_u ** 2 + preri_v ** 2
    preri_R = (np.sqrt(speed_sq) * preri_delta_x3d * 1000) / 1.5e-5
    preri_R[np.isinf(preri_R)] = np.nan
    # Compute the average Reynolds number
    preri_R_avg = np.nanmean(preri_R, axis=0)
    # Compute Reynolds number
    speed_sq = ri_u ** 2 + ri_v ** 2
    ri_R = (np.sqrt(speed_sq) * ri_delta_x3d * 1000) / 1.5e-5
    ri_R[np.isinf(ri_R)] = np.nan
    # Compute the average Reynolds number
    ri_R_avg = np.nanmean(ri_R, axis=0)
    # Compute Reynolds number
    speed_sq = postri_u ** 2 + postri_v ** 2
    postri_R = (np.sqrt(speed_sq) * postri_delta_x3d * 1000) / 1.5e-5
    postri_R[np.isinf(postri_R)] = np.nan
    # Compute the average Reynolds number
    postri_R_avg = np.nanmean(postri_R, axis=0)

    # Plot average Reynolds number
    plt.figure(figsize=(10, 6))
    plt.contourf(preri_R_avg, cmap='nipy_spectral', vmin=0, vmax=1e10,
                 levels=np.linspace(0, 1e10, 20),
                 extent=(np.min(x3d), np.max(x3d), np.min(z3d), np.max(z3d)))
    plt.ylabel('Height (km)')
    plt.xlabel(r'Distance from TC Center (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.title('Average Reynolds Number Before RI')
    plt.colorbar(label='Reynolds Number (unitless)', ticks=np.arange(0, 1.000000001e10, 1e9))
    plt.savefig(os.path.join(base_path, 'preri_rey_avg.png'))
    plt.close()
    # Plot average Reynolds number
    plt.figure(figsize=(10, 6))
    plt.contourf(ri_R_avg, cmap='nipy_spectral', vmin=0, vmax=1e10,
                 levels=np.linspace(0, 1e10, 20),
                 extent=(np.min(x3d), np.max(x3d), np.min(z3d), np.max(z3d)))
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.xlabel(r'Distance from TC Center (km)')
    plt.title('Average Reynolds Number During RI')
    plt.colorbar(label='Reynolds Number (unitless)', ticks=np.arange(0, 1.000000001e10, 1e9))
    plt.savefig(os.path.join(base_path, 'ri_rey_avg.png'))
    plt.close()
    # Plot average Reynolds number
    plt.figure(figsize=(10, 6))
    plt.contourf(postri_R_avg, cmap='nipy_spectral', vmin=0, vmax=1e10,
                 levels=np.linspace(0, 1e10, 20),
                 extent=(np.min(x3d), np.max(x3d), np.min(z3d), np.max(z3d)))
    plt.ylabel('Height (km)')
    plt.xlabel(r'Distance from TC Center (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.title('Average Reynolds Number After RI')
    plt.colorbar(label='Reynolds Number (unitless)', ticks=np.arange(0, 1.000000001e10, 1e9))
    plt.savefig(os.path.join(base_path, 'postri_rey_avg.png'))
    plt.close()
