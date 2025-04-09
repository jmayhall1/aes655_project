# coding=utf-8
"""
Last Edited: 04/09/2025
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np

if __name__ == "__main__":
    # Define file paths
    base_path = '//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/'
    output_path = os.path.join(base_path, 'rey_hourly')
    os.makedirs(output_path, exist_ok=True)
    base_path = '//uahdata/rstor/aes655_project/'

    # Load NetCDF data efficiently
    ncfile = netCDF4.Dataset(os.path.join(base_path, 'cm1out_azimavg_s.nc'))
    z = np.asarray(ncfile.variables['lev'])
    x = np.asarray(ncfile.variables['lon'])
    u = np.asarray(ncfile.variables['u'])[:, :, 0, :]
    v = np.asarray(ncfile.variables['v'])[:, :, 0, :]
    ncfile.close()  # Close file after loading data

    # Compute spatial step differences
    delta_x = np.diff(x, prepend=x[0])

    bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
    top = np.searchsorted(z, 20) + 1

    # Reshape arrays efficiently
    original_zdim, xdim = z.shape[0], x.shape[0]
    z3d = np.tile(z[:, np.newaxis], (1, xdim))[np.newaxis, :, :]
    z3d = np.repeat(z3d, u.shape[0], axis=0)[:, bottom: top, :]

    x3d = np.tile(x[np.newaxis, :], (original_zdim, 1))
    x3d = np.repeat(x3d[np.newaxis, :, :], u.shape[0], axis=0)[:, bottom: top, :]

    delta_x3d = np.tile(delta_x[np.newaxis, :], (original_zdim, 1))
    delta_x3d = np.repeat(delta_x3d[np.newaxis, :, :], u.shape[0], axis=0)[:, bottom: top, :]

    # Slice u and v
    u, v = u[:, bottom: top, :], v[:, bottom: top, :]

    # Compute Reynolds number
    speed_sq = u ** 2 + v ** 2
    R = (np.sqrt(speed_sq) * delta_x3d * 1000) / 1.5e-5
    R[np.isinf(R)] = np.nan

    base_path = '//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Average_Data/Reynolds/'
    # Compute the average Reynolds number
    R_avg = np.nanmean(R, axis=0)

    # Plot average Reynolds number
    plt.figure(figsize=(10, 6))
    plt.contourf(R_avg, cmap='nipy_spectral', vmin=0, vmax=1e10,
                 levels=np.linspace(0, 1e10, 20), extent=(np.min(x3d), np.max(x3d), np.min(z3d), np.max(z3d)))
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.ylabel('Height (km)')
    plt.xlabel(r'Distance from TC Center (km)')
    plt.title('Average Reynolds Number')
    plt.colorbar(label='Reynolds Number (unitless)', ticks=np.arange(0, 1.000000001e10, 1e9))
    plt.savefig(os.path.join(base_path, 'rey_avg.png'), dpi=300)
    plt.close()

    # Save individual timestep images efficiently
    for i, Ri in enumerate(R):
        plt.figure(figsize=(10, 6))
        plt.contourf(Ri, cmap='nipy_spectral', vmin=0, vmax=1e10,
                     levels=np.linspace(0, 1e10, 20),
                     extent=(np.min(x3d), np.max(x3d), np.min(z3d), np.max(z3d)))
        plt.ylabel('Height (km)')
        plt.yticks(ticks=np.arange(10, 21, 2))
        plt.xticks(ticks=np.arange(0, 301, 50))
        plt.xlabel(r'Distance from TC Center (km)')
        plt.title(f'Reynolds Number at Hour {i}')
        plt.colorbar(label='Reynolds Number (unitless)', ticks=np.arange(0, 1.000000001e10, 1e9))
        plt.savefig(os.path.join(output_path, f'rey_{i}.png'), dpi=300)
        plt.close()
