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
    base_path = '//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/'
    output_path = os.path.join(base_path, 'rey_hourly')
    os.makedirs(output_path, exist_ok=True)
    base_path = '//uahdata/rstor/aes655_project/'

    # Load NetCDF data efficiently
    ncfile = netCDF4.Dataset(os.path.join(base_path, 'cm1out_azimavg_s.nc'))
    z = np.array(ncfile.variables['lev'])
    x = np.array(ncfile.variables['lon'])
    u = np.array(ncfile.variables['u'])[:, :, 0, :]
    v = np.array(ncfile.variables['v'])[:, :, 0, :]
    ncfile.close()  # Close file after loading data

    # Compute spatial step differences
    delta_x = np.diff(x, prepend=x[0])

    # Find the bottom index where z > 10
    bottom = np.searchsorted(z, 10)

    # Reshape arrays efficiently
    original_zdim, xdim = z.shape[0], x.shape[0]
    z3d = np.tile(z[:, np.newaxis], (1, xdim))[np.newaxis, :, :]
    z3d = np.repeat(z3d, u.shape[0], axis=0)[:, bottom:, :]

    x3d = np.tile(x[np.newaxis, :], (original_zdim, 1))
    x3d = np.repeat(x3d[np.newaxis, :, :], u.shape[0], axis=0)[:, bottom:, :]

    delta_x3d = np.tile(delta_x[np.newaxis, :], (original_zdim, 1))
    delta_x3d = np.repeat(delta_x3d[np.newaxis, :, :], u.shape[0], axis=0)[:, bottom:, :]

    # Slice u and v
    u, v = u[:, bottom:, :], v[:, bottom:, :]

    # Compute Reynolds number
    speed_sq = u ** 2 + v ** 2
    R = (np.sqrt(speed_sq) * delta_x3d * 1000) / 1.5e-5
    R[np.isinf(R)] = np.nan

    base_path = '//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Average_Data/Reynolds/'
    # Compute the average Reynolds number
    R_avg = np.nanmean(R, axis=0)

    # Plot average Reynolds number
    plt.figure(figsize=(8, 6))
    plt.imshow(R_avg, aspect='auto', cmap='nipy_spectral', vmin=0,
               extent=(np.min(x3d), np.max(x3d), np.max(z3d), np.min(z3d)))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel('Range (km)')
    plt.title('Average Reynolds Number')
    plt.colorbar(label='Reynolds Number (unitless)')
    plt.savefig(os.path.join(base_path, 'rey_avg.png'))
    plt.close()

    # Save individual timestep images efficiently
    for i, Ri in enumerate(R):
        plt.figure(figsize=(8, 6))
        plt.imshow(Ri, aspect='auto', cmap='nipy_spectral', vmin=0,
                   extent=(np.min(x3d), np.max(x3d), np.max(z3d) / 1000, np.min(z3d) / 1000))
        plt.gca().invert_yaxis()
        plt.ylabel('Height (km)')
        plt.xlabel('Range (km)')
        plt.title(f'Reynolds Number at Timestep {i}')
        plt.colorbar(label='Reynolds Number (unitless)')
        plt.savefig(os.path.join(output_path, f'rey_{i}.png'))
        plt.close()
