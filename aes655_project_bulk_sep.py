# coding=utf-8
"""
Last Edited: 04/09/2025
@author: John Mark Mayhall
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import netCDF4
import numpy as np


# Function to save occurrence heatmaps
def save_heatmap(data: np.array, title: str, filename: str, colorbar_label: str, cticks: np.array,
                 levels: np.array) -> None:
    """
    Function to plot occurrence heatmaps.
    :param data: Data to be used for plotting
    :param title: Title of the created plot
    :param filename: Saved plot's filename
    :param colorbar_label: Colorbar Label
    :param cticks: Colorbar ticks
    :param levels: Plot color levels
    :return: Nothing
    """
    plt.figure(figsize=(10, 6))
    plt.contourf(data, cmap='turbo', vmin=0,
                 extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000), levels=levels)
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.xlabel(r'Distance from TC Center (km)')
    plt.title(title)
    plt.colorbar(label=colorbar_label, ticks=cticks)
    plt.savefig(os.path.join(base_path, filename), dpi=300)
    plt.close()


def save_heatmap_double(data1: np.array, data2: np.array, title: str, filename: str, colorbar_label: str,
                        cticks: np.array, levels: np.array) -> None:
    """
    Function to plot occurrence heatmaps.
    :param data1: Data to be used for plotting
    :param data2: Data to be used for plotting
    :param title: Title of the created plot
    :param filename: Saved plot's filename
    :param colorbar_label: Colorbar Label
    :param cticks: Colorbar ticks
    :param levels: Plot color levels
    :return: Nothing
    """
    fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
    im1 = ax.contourf(data2, cmap='terrain', vmin=0, vmax=60,
                      extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000),
                      alpha=0.75, levels=np.linspace(0, 60, 12))
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 60), cmap='terrain'),
                 ax=ax, orientation='vertical', label=colorbar_label, ticks=np.linspace(0, 60, 12))
    im2 = ax.contour(data1, cmap='turbo', vmin=0, vmax=1,
                     extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000), levels=levels, linewidths=3)
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap='turbo'),
                 ax=ax, orientation='vertical', label=(r'# of $R_b$ < 1 Occurrences' + '\n Divided by # of Timesteps'),
                 ticks=levels)
    ax.set_ylabel('Height (km)')
    ax.set_yticks(ticks=np.arange(10, 21, 2))
    ax.set_xticks(ticks=np.arange(0, 301, 50))
    ax.set_xlabel(r'Distance from TC Center (km)')
    ax.set_title(title)
    plt.savefig(os.path.join(base_path, filename), dpi=300)
    plt.close()


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
    post_ri_time = np.where(mtime == 76)[0][0]
    preri_u = np.asarray(ncfile.variables['u'])[:pre_ri_time, :, 0, :]
    preri_v = np.asarray(ncfile.variables['v'])[:pre_ri_time, :, 0, :]
    preri_thv = np.asarray(ncfile.variables['thv'])[:pre_ri_time, :, 0, :]
    ri_u = np.asarray(ncfile.variables['u'])[pre_ri_time:post_ri_time, :, 0, :]
    ri_v = np.asarray(ncfile.variables['v'])[pre_ri_time:post_ri_time, :, 0, :]
    ri_thv = np.asarray(ncfile.variables['thv'])[pre_ri_time:post_ri_time, :, 0, :]
    postri_u = np.asarray(ncfile.variables['u'])[post_ri_time:, :, 0, :]
    postri_v = np.asarray(ncfile.variables['v'])[post_ri_time:, :, 0, :]
    postri_thv = np.asarray(ncfile.variables['thv'])[post_ri_time:, :, 0, :]
    ncfile.close()  # Close file after loading data

    # Compute vertical differences
    delta_z = np.diff(z, prepend=z[0])

    # Compute vertical derivatives using MetPy
    preri_delta_u = mpcalc.first_derivative(preri_u, axis=1, x=z)
    preri_delta_v = mpcalc.first_derivative(preri_v, axis=1, x=z)
    preri_delta_thv = mpcalc.first_derivative(preri_thv, axis=1, x=z)
    # Compute vertical derivatives using MetPy
    ri_delta_u = mpcalc.first_derivative(ri_u, axis=1, x=z)
    ri_delta_v = mpcalc.first_derivative(ri_v, axis=1, x=z)
    ri_delta_thv = mpcalc.first_derivative(ri_thv, axis=1, x=z)
    # Compute vertical derivatives using MetPy
    postri_delta_u = mpcalc.first_derivative(postri_u, axis=1, x=z)
    postri_delta_v = mpcalc.first_derivative(postri_v, axis=1, x=z)
    postri_delta_thv = mpcalc.first_derivative(postri_thv, axis=1, x=z)

    bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
    top = np.searchsorted(z, 20) + 1

    # Trim arrays to remove unnecessary levels
    z = np.tile(z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    z = np.repeat(z, preri_u.shape[0], axis=0)[:, bottom: top, :] * 1000

    preri_delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    preri_delta_z = np.repeat(preri_delta_z, preri_u.shape[0], axis=0)[:, bottom: top, :] * 1000
    ri_delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    ri_delta_z = np.repeat(ri_delta_z, ri_u.shape[0], axis=0)[:, bottom: top, :] * 1000
    postri_delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    postri_delta_z = np.repeat(postri_delta_z, postri_u.shape[0], axis=0)[:, bottom: top, :] * 1000

    preri_u, preri_v, preri_thv = preri_u[:, bottom: top, :], preri_v[:, bottom: top, :], preri_thv[:, bottom: top, :]
    preri_delta_u, preri_delta_v, preri_delta_thv = (preri_delta_u[:, bottom: top, :], preri_delta_v[:, bottom: top, :],
                                                     preri_delta_thv[:, bottom: top, :])
    ri_u, ri_v, ri_thv = ri_u[:, bottom: top, :], ri_v[:, bottom: top, :], ri_thv[:, bottom: top, :]
    ri_delta_u, ri_delta_v, ri_delta_thv = ri_delta_u[:, bottom: top, :], ri_delta_v[:, bottom: top, :], ri_delta_thv[:,
                                                                                                         bottom: top, :]
    postri_u, postri_v, postri_thv = (postri_u[:, bottom: top, :], postri_v[:, bottom: top, :],
                                      postri_thv[:, bottom: top, :])
    postri_delta_u, postri_delta_v, postri_delta_thv = (postri_delta_u[:, bottom: top, :],
                                                        postri_delta_v[:, bottom: top, :],
                                                        postri_delta_thv[:, bottom: top, :])

    # Compute Richardson number
    preri_R = (9.81 * preri_delta_z * preri_delta_thv) / (preri_thv * (preri_delta_u ** 2 + preri_delta_v ** 2))
    preri_R[np.isinf(preri_R)] = np.nan
    # Compute the average Richardson number
    preri_R_avg = np.nanmean(preri_R, axis=0)
    # Compute Richardson number
    ri_R = (9.81 * ri_delta_z * ri_delta_thv) / (ri_thv * (ri_delta_u ** 2 + ri_delta_v ** 2))
    ri_R[np.isinf(ri_R)] = np.nan
    # Compute the average Richardson number
    ri_R_avg = np.nanmean(ri_R, axis=0)
    # Compute Richardson number
    postri_R = (9.81 * postri_delta_z * postri_delta_thv) / (postri_thv * (postri_delta_u ** 2 + postri_delta_v ** 2))
    postri_R[np.isinf(postri_R)] = np.nan
    # Compute the average Richardson number
    postri_R_avg = np.nanmean(postri_R, axis=0)

    preri_wind_avg = np.nanmean(np.sqrt(preri_u ** 2 + preri_v ** 2), axis=0)
    ri_wind_avg = np.nanmean(np.sqrt(ri_u ** 2 + ri_v ** 2), axis=0)
    postri_wind_avg = np.nanmean(np.sqrt(postri_u ** 2 + postri_v ** 2), axis=0)

    # Save the average Richardson number plot
    preri_R_avg[preri_R_avg > 10] = 10
    plt.figure(figsize=(10, 6))
    plt.contourf(preri_R_avg, vmin=-1, vmax=10, cmap='turbo', levels=np.arange(0, 11, 1),
                 extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000))
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.xlabel('Range (km)')
    plt.title('Average Richardson Number Before RI')
    plt.colorbar(label=r'$R_b$ (unitless)', ticks=np.arange(0, 11, 2))
    plt.savefig(os.path.join(base_path, 'preri_Rb_avg.png'), dpi=300)
    plt.close()
    # Save the average Richardson number plot
    ri_R_avg[ri_R_avg > 10] = 10
    plt.figure(figsize=(10, 6))
    plt.contourf(ri_R_avg, cmap='turbo', vmin=-1, vmax=10, levels=np.arange(0, 11, 1),
                 extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000))
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.xlabel('Range (km)')
    plt.title('Average Richardson Number During RI')
    plt.colorbar(label=r'$R_b$ (unitless)', ticks=np.arange(0, 11, 2))
    plt.savefig(os.path.join(base_path, 'ri_Rb_avg.png'), dpi=300)
    plt.close()
    # Save the average Richardson number plot
    postri_R_avg[postri_R_avg > 10] = 10
    plt.figure(figsize=(10, 6))
    plt.contourf(postri_R_avg, cmap='turbo', vmin=-1, vmax=10, levels=np.arange(0, 11, 1),
                 extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000))
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.ylabel('Height (km)')
    plt.xlabel('Range (km)')
    plt.title('Average Richardson Number After RI')
    plt.colorbar(label=r'$R_b$ (unitless)', ticks=np.arange(0, 11, 2))
    plt.savefig(os.path.join(base_path, 'postri_Rb_avg.png'), dpi=300)
    plt.close()

    # Initialize tracking arrays
    num_timesteps, zdim, xdim = preri_R.shape
    minimum_loc = np.zeros((zdim, xdim), dtype=int)
    loc025 = np.zeros((zdim, xdim), dtype=int)
    loc1 = np.zeros((zdim, xdim), dtype=int)
    # Loop through timesteps for visualization and tracking
    for i in range(num_timesteps):
        Ri = preri_R[i, :, :]

        # Identify critical regions
        minimum_loc += (Ri == np.nanmin(Ri)).astype(int)
        loc025 += (Ri < 0.25).astype(int)
        loc1 += (Ri < 1).astype(int)
    minimum_loc = np.divide(minimum_loc, num_timesteps)
    loc025 = np.divide(loc025, num_timesteps)
    loc1 = np.divide(loc1, num_timesteps)
    # Save occurrence heatmaps
    save_heatmap(minimum_loc, 'Location of Minimum $R_b$ Before RI \n Normalized By Number of Timesteps',
                 'preri_minRb_loc.png',
                 colorbar_label=(r'# of Minimum $R_b$ Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 0.2, 10), levels=np.linspace(0, 0.2, 10))
    save_heatmap(loc025, 'Location of $R_b$ < 0.25 Before RI \n Normalized By Number of Timesteps',
                 'preri_025Rb_loc.png',
                 colorbar_label=(r'# of $R_b$ < 0.25 Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 1, 10), levels=np.linspace(0, 1, 20))
    save_heatmap_double(data1=loc025, data2=preri_wind_avg,
                        title='Location of $R_b$ < 0.25 Before RI \n Normalized By Number of'
                              'Timesteps Contoured over Velocity Before RI',
                        filename='preri_025Rb_wind_loc.png',
                        colorbar_label=(r'Velocity ($\frac{m}{s}$)'),
                        cticks=np.linspace(0, 1, 10), levels=np.linspace(0, 1, 20))
    save_heatmap(loc1, 'Location of $R_b$ < 1 Before RI \n Normalized By Number of Timesteps',
                 'preri_1Rb_loc.png',
                 colorbar_label=(r'# of $R_b$ < 1 Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 1.1, 10), levels=np.linspace(0, 1.1, 20))

    # Initialize tracking arrays
    num_timesteps, zdim, xdim = ri_R.shape
    minimum_loc = np.zeros((zdim, xdim), dtype=int)
    loc025 = np.zeros((zdim, xdim), dtype=int)
    loc1 = np.zeros((zdim, xdim), dtype=int)
    # Loop through timesteps for visualization and tracking
    for i in range(num_timesteps):
        Ri = ri_R[i, :, :]

        # Identify critical regions
        minimum_loc += (Ri == np.nanmin(Ri)).astype(int)
        loc025 += (Ri < 0.25).astype(int)
        loc1 += (Ri < 1).astype(int)
    minimum_loc = np.divide(minimum_loc, num_timesteps)
    loc025 = np.divide(loc025, num_timesteps)
    loc1 = np.divide(loc1, num_timesteps)
    # Save occurrence heatmaps
    save_heatmap(minimum_loc, 'Location of Minimum $R_b$ During RI \n Normalized By Number of Timesteps',
                 'ri_minRb_loc.png',
                 colorbar_label=(r'# of Minimum $R_b$ Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 0.2, 10), levels=np.linspace(0, 0.2, 10))
    save_heatmap(loc025, 'Location of $R_b$ < 0.25 During RI \n Normalized By Number of Timesteps',
                 'ri_025Rb_loc.png',
                 colorbar_label=(r'# of $R_b$ < 0.25 Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 1, 10), levels=np.linspace(0, 1, 20))
    save_heatmap_double(data1=loc025, data2=ri_wind_avg,
                        title='Location of $R_b$ < 0.25 During RI \n Normalized By Number of'
                              'Timesteps Contoured over During RI',
                        filename='ri_025Rb_wind_loc.png',
                        colorbar_label=(r'Velocity ($\frac{m}{s}$)'),
                        cticks=np.linspace(0, 1, 10), levels=np.linspace(0, 1, 20))
    save_heatmap(loc1, 'Location of $R_b$ < 1 During RI \n Normalized By Number of Timesteps',
                 'ri_1Rb_loc.png',
                 colorbar_label=(r'# of $R_b$ < 1 Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 1.1, 10), levels=np.linspace(0, 1.1, 20))

    # Initialize tracking arrays
    num_timesteps, zdim, xdim = postri_R.shape
    minimum_loc = np.zeros((zdim, xdim), dtype=int)
    loc025 = np.zeros((zdim, xdim), dtype=int)
    loc1 = np.zeros((zdim, xdim), dtype=int)
    # Loop through timesteps for visualization and tracking
    for i in range(num_timesteps):
        Ri = postri_R[i, :, :]

        # Identify critical regions
        minimum_loc += (Ri == np.nanmin(Ri)).astype(int)
        loc025 += (Ri < 0.25).astype(int)
        loc1 += (Ri < 1).astype(int)
    minimum_loc = np.divide(minimum_loc, num_timesteps)
    loc025 = np.divide(loc025, num_timesteps)
    loc1 = np.divide(loc1, num_timesteps)
    # Save occurrence heatmaps
    save_heatmap(minimum_loc, 'Location of Minimum $R_b$ After RI \n Normalized By Number of Timesteps',
                 'postri_minRb_loc.png',
                 colorbar_label=(r'# of Minimum $R_b$ Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 0.2, 10), levels=np.linspace(0, 0.2, 10))
    save_heatmap(loc025, 'Location of $R_b$ < 0.25 After RI \n Normalized By Number of Timesteps',
                 'postri_025Rb_loc.png',
                 colorbar_label=(r'# of $R_b$ < 0.25 Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 1, 10), levels=np.linspace(0, 1, 20))
    save_heatmap_double(data1=loc025, data2=postri_wind_avg,
                        title='Location of $R_b$ < 0.25 After RI \n Normalized By Number of'
                              'Timesteps Contoured over Velocity after RI',
                        filename='postri_025Rb_wind_loc.png',
                        colorbar_label=(r'Velocity ($\frac{m}{s}$)'),
                        cticks=np.linspace(0, 1, 10), levels=np.linspace(0, 1, 20))
    save_heatmap(loc1, 'Location of $R_b$ < 1 After RI \n Normalized By Number of Timesteps',
                 'postri_1Rb_loc.png',
                 colorbar_label=(r'# of $R_b$ < 1 Occurrences' + '\n Divided by # of Timesteps'),
                 cticks=np.linspace(0, 1.1, 10), levels=np.linspace(0, 1.1, 20))
