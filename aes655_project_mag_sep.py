# coding=utf-8
"""
Last Edited: 04/09/2025
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np


# Function to generate plots
def plot_and_save(data: np.array, title: str, filename: str) -> None:
    """
    Function for plotting data
    :param data: Array to be plotted
    :param title: Title of plot.
    :param filename: Filename of the plot that is being saved.
    :return: Nothing
    """
    plt.figure(figsize=(10, 6))
    corrected_data = np.copy(data)
    corrected_data[corrected_data > 110] = 110
    plt.contourf(corrected_data, cmap='turbo', vmin=0, vmax=110, levels=np.arange(0, 111, 10),
                 extent=(x.min(), x.max(), z.min(), z.max()))
    plt.ylabel('Height (km)')
    plt.xlabel(r'Distance from TC Center (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.title(title)
    plt.colorbar(label='%', ticks=np.arange(0, 101, 20))
    plt.savefig(f'//uahdata/rstor/aes655_project/sep_by_intensity_phase/{filename}.png', dpi=300)
    plt.close()


# Function for safe division to handle inf values
def safe_divide(numerator: np.array, denominator: np.array) -> np.array:
    """
    Function for dividing arrays to find percents.
    :param numerator: Array that is being compared.
    :param denominator: Array that is being compared to.
    :return: Array of Percentages
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.abs(numerator / denominator) * 100
        result[np.isinf(result)] = np.nan
    return result


if __name__ == "__main__":
    # The file path
    path = '//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc'
    # Load netCDF dataset once
    dataset = netCDF4.Dataset(path)

    # Read variables
    z = np.asarray(dataset.variables['lev'])
    x = np.asarray(dataset.variables['lon'])
    bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
    top = np.searchsorted(z, 20) + 1
    z = z[bottom: top]
    mtime = (np.asarray(dataset.variables['mtime']) / 3600)[:, 0, 0]
    pre_ri_time = np.where(mtime == 10)[0][0]
    post_ri_time = np.where(mtime == 76)[0][0]

    # Extract data slices only once
    data_shear = np.asarray(dataset.variables['qshear'])[:pre_ri_time, bottom: top, 0, :]
    data_buoy = np.asarray(dataset.variables['qbuoy'])[:pre_ri_time, bottom: top, 0, :]
    data_diss = np.asarray(dataset.variables['qdiss'])[:pre_ri_time, bottom: top, 0, :]
    # Compute normalized values
    buoy_shear = safe_divide(data_buoy, data_shear)
    diss_shear = safe_divide(data_diss, data_shear)
    buoydiss_shear = safe_divide(data_buoy + data_diss, data_shear)
    # Compute averages
    buoy_avg_shear = np.nanmean(buoy_shear, axis=0)
    diss_avg_shear = np.nanmean(diss_shear, axis=0)
    buoydiss_avg_shear = np.nanmean(buoydiss_shear, axis=0)
    # Save average plots
    plot_and_save(buoy_avg_shear, 'Average Buoyancy Production % Before RI', 'preri_buoy_avg_percent')
    plot_and_save(diss_avg_shear, 'Average Dissipation % Before RI', 'preri_diss_avg_percent')
    plot_and_save(buoydiss_avg_shear, 'Average Buoyancy Production and Dissipation % Before RI',
                  'preri_buoydiss_avg_percent')

    # Extract data slices only once
    data_shear = np.asarray(dataset.variables['qshear'])[pre_ri_time:post_ri_time, bottom: top, 0, :]
    data_buoy = np.asarray(dataset.variables['qbuoy'])[pre_ri_time:post_ri_time, bottom: top, 0, :]
    data_diss = np.asarray(dataset.variables['qdiss'])[pre_ri_time:post_ri_time, bottom: top, 0, :]
    # Compute normalized values
    buoy_shear = safe_divide(data_buoy, data_shear)
    diss_shear = safe_divide(data_diss, data_shear)
    buoydiss_shear = safe_divide(data_buoy + data_diss, data_shear)
    # Compute averages
    buoy_avg_shear = np.nanmean(buoy_shear, axis=0)
    diss_avg_shear = np.nanmean(diss_shear, axis=0)
    buoydiss_avg_shear = np.nanmean(buoydiss_shear, axis=0)
    # Save average plots
    plot_and_save(buoy_avg_shear, 'Average Buoyancy Production % During RI', 'ri_buoy_avg_percent')
    plot_and_save(diss_avg_shear, 'Average Dissipation % During RI', 'ri_diss_avg_percent')
    plot_and_save(buoydiss_avg_shear, 'Average Buoyancy Production and Dissipation % During RI',
                  'ri_buoydiss_avg_percent')

    # Extract data slices only once
    data_shear = np.asarray(dataset.variables['qshear'])[post_ri_time:, bottom: top, 0, :]
    data_buoy = np.asarray(dataset.variables['qbuoy'])[post_ri_time:, bottom: top, 0, :]
    data_diss = np.asarray(dataset.variables['qdiss'])[post_ri_time:, bottom: top, 0, :]
    # Compute normalized values
    buoy_shear = safe_divide(data_buoy, data_shear)
    diss_shear = safe_divide(data_diss, data_shear)
    buoydiss_shear = safe_divide(data_buoy + data_diss, data_shear)
    # Compute averages
    buoy_avg_shear = np.nanmean(buoy_shear, axis=0)
    diss_avg_shear = np.nanmean(diss_shear, axis=0)
    buoydiss_avg_shear = np.nanmean(buoydiss_shear, axis=0)
    # Save average plots
    plot_and_save(buoy_avg_shear, 'Average Buoyancy Production % After RI', 'postri_buoy_avg_percent')
    plot_and_save(diss_avg_shear, 'Average Dissipation % After RI', 'postri_diss_avg_percent')
    plot_and_save(buoydiss_avg_shear, 'Average Buoyancy Production and Dissipation % After RI',
                  'postri_buoydiss_avg_percent')
