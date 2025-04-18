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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm


# Function to plot data
def plot_and_save(data_func: np.array, title: str, filename: str, extent: tuple, xticks: np.array,
                  cmap: str, cticks: np.array, levels: np.array, cbar_label='TKE per Second ($m^2/s^3$)') -> None:
    """
    Function for plotting arrays.
    :param norm: Colorbar normalization
    :param cticks: Colorbar ticks
    :param xticks: X-axis ticks
    :param levels: Colorbar levels
    :param data_func: Array to be plotted.
    :param title: Title of plot.
    :param filename: Filename of plot being saved.
    :param cmap: Colormap of plot.
    :param cbar_label: Label of the colorbar.
    :param extent: The extent.
    :return: Nothing.
    """
    plt.figure(figsize=(10, 6))
    plt.contourf(data_func, cmap=cmap, extent=extent, levels=levels)
    plt.xlabel(r'Distance from TC Center (km)')
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(xticks)
    plt.title(title)
    plt.colorbar(label=cbar_label, ticks=cticks)
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()


def compute_dominant_term(*terms):
    """
    Calculate dominate production term
    :param terms: Terms needed for dominant term calculation.
    :return: Array with integers representing the max term.
    """
    abs_terms = [np.abs(term) for term in terms]
    dominant = np.argmax(abs_terms, axis=0) + 1  # Add 1 to match colormap labels
    return dominant

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
    im2 = ax.contour(data1, cmap='turbo', vmin=0, vmax=0.002,
                     extent=(np.min(x), np.max(x), np.min(z) / 1000, np.max(z) / 1000), levels=levels, linewidths=3)
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 0.002), cmap='turbo'),
                 ax=ax, orientation='vertical', label=(r'$N^2 (\frac{1}{s})$'),
                 ticks=levels)
    ax.set_ylabel('Height (km)')
    ax.set_yticks(ticks=np.arange(10, 21, 2))
    ax.set_xticks(ticks=np.arange(0, 301, 50))
    ax.set_xlabel(r'Distance from TC Center (km)')
    ax.set_title(title)
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close()


if __name__ == "__main__":
    # Define the file path
    path = "//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc"
    output_dir = "//uahdata/rstor/aes655_project/sep_by_intensity_phase/"

    # Load dataset
    with netCDF4.Dataset(path) as ds:
        z = np.asarray(ds.variables['lev'])
        x = np.asarray(ds.variables['lon'])
        mtime = (np.asarray(ds.variables['mtime']) / 3600)[:, 0, 0]
        pre_ri_time = np.where(mtime == 10)[0][0]
        post_ri_time = np.where(mtime == 76)[0][0]
        bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
        top = np.searchsorted(z, 20) + 1
        z = z[bottom: top]
        extent = (np.min(x), np.max(x), np.min(z), np.max(z))

        # Extract and slice data
        variables = ['rho', 'thv']
        preri_data = {var: np.asarray(ds.variables[var])[:pre_ri_time, bottom: top, 0, :] for var in variables}
        ri_data = {var: np.asarray(ds.variables[var])[pre_ri_time:post_ri_time, bottom: top, 0, :] for var in variables}
        postri_data = {var: np.asarray(ds.variables[var])[post_ri_time:, bottom: top, 0, :] for var in variables}

    # Compute average values
    preri_avg_data = {key: np.nanmean(val, axis=0) for key, val in preri_data.items()}
    ri_avg_data = {key: np.nanmean(val, axis=0) for key, val in ri_data.items()}
    postri_avg_data = {key: np.nanmean(val, axis=0) for key, val in postri_data.items()}

    path = "//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc"
    output_dir = "//uahdata/rstor/aes655_project/sep_by_intensity_phase/"
    base_path = '//uahdata/rstor/aes655_project/'
    ncfile = netCDF4.Dataset(os.path.join(base_path, 'cm1out_azimavg_s.nc'))
    z = np.asarray(ncfile.variables['lev'])
    x = np.asarray(ncfile.variables['lon'])

    preri_u = np.asarray(ncfile.variables['u'])[:pre_ri_time, :, 0, :]
    preri_v = np.asarray(ncfile.variables['v'])[:pre_ri_time, :, 0, :]
    preri_u, preri_v = preri_u[:, bottom: top, :], preri_v[:, bottom: top, :]
    ri_u = np.asarray(ncfile.variables['u'])[pre_ri_time:post_ri_time, :, 0, :]
    ri_v = np.asarray(ncfile.variables['v'])[pre_ri_time:post_ri_time, :, 0, :]
    ri_u, ri_v = ri_u[:, bottom: top, :], ri_v[:, bottom: top, :]
    postri_u = np.asarray(ncfile.variables['u'])[post_ri_time:, :, 0, :]
    postri_v = np.asarray(ncfile.variables['v'])[post_ri_time:, :, 0, :]
    postri_u, postri_v = postri_u[:, bottom: top, :], postri_v[:, bottom: top, :]
    preri_wind_avg = np.nanmean(np.sqrt(preri_u ** 2 + preri_v ** 2), axis=0)
    ri_wind_avg = np.nanmean(np.sqrt(ri_u ** 2 + ri_v ** 2), axis=0)
    postri_wind_avg = np.nanmean(np.sqrt(postri_u ** 2 + postri_v ** 2), axis=0)

    mtime = (np.asarray(ncfile.variables['mtime']) / 3600)[:, 0, 0]
    pre_ri_time = np.where(mtime == 10)[0][0]
    post_ri_time = np.where(mtime == 76)[0][0]
    preri_rho = np.asarray(ncfile.variables['rho'])[:pre_ri_time, :, 0, :]
    preri_thv = np.asarray(ncfile.variables['thv'])[:pre_ri_time, :, 0, :]
    ri_rho = np.asarray(ncfile.variables['rho'])[pre_ri_time:post_ri_time, :, 0, :]
    ri_thv = np.asarray(ncfile.variables['thv'])[pre_ri_time:post_ri_time, :, 0, :]
    postri_rho = np.asarray(ncfile.variables['rho'])[post_ri_time:, :, 0, :]
    postri_thv = np.asarray(ncfile.variables['thv'])[post_ri_time:, :, 0, :]
    ncfile.close()  # Close file after loading data

    # Compute vertical differences
    delta_z = np.diff(z, prepend=z[0])

    # Compute vertical derivatives using MetPy
    preri_delta_thv = mpcalc.first_derivative(preri_thv, axis=1, x=z)
    ri_delta_thv = mpcalc.first_derivative(ri_thv, axis=1, x=z)
    postri_delta_thv = mpcalc.first_derivative(postri_thv, axis=1, x=z)

    bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
    top = np.searchsorted(z, 20) + 1

    # Trim arrays to remove unnecessary levels
    z = np.tile(z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    z = np.repeat(z, preri_rho.shape[0], axis=0)[:, bottom: top, :] * 1000

    preri_delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    preri_delta_z = np.repeat(preri_delta_z, preri_rho.shape[0], axis=0)[:, bottom: top, :] * 1000
    ri_delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    ri_delta_z = np.repeat(ri_delta_z, ri_rho.shape[0], axis=0)[:, bottom: top, :] * 1000
    postri_delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
    postri_delta_z = np.repeat(postri_delta_z, postri_rho.shape[0], axis=0)[:, bottom: top, :] * 1000

    preri_rho, preri_thv = preri_rho[:, bottom: top, :], preri_thv[:, bottom: top, :]
    preri_delta_thv = preri_delta_thv[:, bottom: top, :]
    ri_rho, ri_thv = ri_rho[:, bottom: top, :], ri_thv[:, bottom: top, :]
    ri_delta_thv = ri_delta_thv[:, bottom: top, :]
    postri_thv = postri_thv[:, bottom: top, :]
    postri_delta_thv = postri_delta_thv[:, bottom: top, :]

    preri_n = np.multiply(np.divide(9.81, preri_thv), np.divide(preri_delta_thv, preri_delta_z))
    preri_n_avg = np.nanmean(preri_n, axis=0)
    ri_n = np.multiply(np.divide(9.81, ri_thv), np.divide(ri_delta_thv, ri_delta_z))
    ri_n_avg = np.nanmean(ri_n, axis=0)
    postri_n = np.multiply(np.divide(9.81, postri_thv), np.divide(postri_delta_thv, postri_delta_z))
    postri_n_avg = np.nanmean(postri_n, axis=0)

    # Plot averaged data
    plot_and_save(preri_avg_data['rho'], 'Average Dry Air Density Before RI', 'preri_rho_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(0, 0.5, 10), cticks=np.linspace(0, 0.5, 10),
                  xticks=np.arange(0, 301, 50), cbar_label=r'$\rho$ $(\frac{kg}{m^3})$')
    plot_and_save(preri_avg_data['thv'], 'Average Virtual Potential Temperature Before RI', 'preri_thv_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(340, 500, 16), cticks=np.linspace(340, 500, 8),
                  xticks=np.arange(0, 301, 50), cbar_label=r'T (K)')
    plot_and_save(preri_n_avg, 'Average Squared Brunt–Väisälä Frequency Before RI', 'preri_n_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(0, 0.002, 10), cticks=np.linspace(0, 0.002, 10),
                  xticks=np.arange(0, 301, 50), cbar_label=r'$N^2 (\frac{1}{s})$')
    save_heatmap_double(data1=preri_n_avg, data2=preri_wind_avg,
                        title='Average Brunt–Väisälä Frequency Contoured over Velocity Before RI',
                        filename='preri_n_wind_loc.png',
                        colorbar_label=(r'Velocity ($\frac{m}{s}$)'),
                        cticks=np.linspace(0, 0.002, 10), levels=np.linspace(0, 0.002, 10))

    # Plot averaged data
    plot_and_save(ri_avg_data['rho'], 'Average Dry Air Density During RI', 'ri_rho_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(0, 0.5, 10), cticks=np.linspace(0, 0.5, 10),
                  xticks=np.arange(0, 301, 50), cbar_label=r'$\rho$ $(\frac{kg}{m^3})$')
    plot_and_save(ri_avg_data['thv'], 'Average Virtual Potential Temperature During RI', 'ri_thv_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(340, 500, 16), cticks=np.arange(340, 500, 8),
                  xticks=np.linspace(0, 301, 50), cbar_label=r'T (K)')
    plot_and_save(ri_n_avg, 'Average Squared Brunt–Väisälä Frequency During RI', 'ri_n_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(0, 0.002, 10), cticks=np.linspace(0, 0.002, 10),
                  xticks=np.arange(0, 301, 50), cbar_label=r'$N^2 (\frac{1}{s})$')
    save_heatmap_double(data1=ri_n_avg, data2=ri_wind_avg,
                        title='Average Brunt–Väisälä Frequency Contoured over Velocity During RI',
                        filename='ri_n_wind_loc.png',
                        colorbar_label=(r'Velocity ($\frac{m}{s}$)'),
                        cticks=np.linspace(0, 0.002, 10), levels=np.linspace(0, 0.002, 10))

    # Plot averaged data
    plot_and_save(postri_avg_data['rho'], 'Average Dry Air Density After RI', 'postri_rho_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(0, 0.5, 10), cticks=np.linspace(0, 0.5, 10),
                  xticks=np.arange(0, 301, 50), cbar_label=r'$\rho$ $(\frac{kg}{m^3})$')
    plot_and_save(postri_avg_data['thv'], 'Average Virtual Potential Temperature After RI', 'postri_thv_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(340, 500, 16), cticks=np.linspace(340, 500, 8),
                  xticks=np.arange(0, 301, 50), cbar_label=r'T (K)')
    plot_and_save(postri_n_avg, 'Average Brunt–Väisälä Frequency After RI', 'postri_n_avg.png',
                  extent, cmap='turbo',
                  levels=np.linspace(0, 0.002, 10), cticks=np.linspace(0, 0.002, 10),
                  xticks=np.arange(0, 301, 50), cbar_label=r'$N^2 (\frac{1}{s})$')
    save_heatmap_double(data1=postri_n_avg, data2=postri_wind_avg,
                        title='Average Brunt–Väisälä Frequency Contoured over Velocity After RI',
                        filename='postri_n_wind_loc.png',
                        colorbar_label=(r'Velocity ($\frac{m}{s}$)'),
                        cticks=np.linspace(0, 0.002, 10), levels=np.linspace(0, 0.002, 10))
