# coding=utf-8
"""
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
    plt.imshow(data, aspect='auto', cmap='rainbow', vmin=0, vmax=110,
               extent=(x.min(), x.max(), z.max(), z.min()))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel('Range (km)')
    plt.title(title)
    plt.colorbar(label='%')
    plt.savefig(f'//uahdata/rstor/aes655_project/sep_by_intensity_phase/{filename}.png')
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
    z = np.array(dataset.variables['lev'])
    x = np.array(dataset.variables['lon'])
    bottom = np.searchsorted(z, 10)  # More efficient than np.where(z > 10)[0][0]
    z = z[bottom:]
    mtime = (np.array(dataset.variables['mtime']) / 3600)[:, 0, 0]
    pre_ri_time = np.where(mtime == 40)[0][0]
    post_ri_time = np.where(mtime == 61)[0][0]

    # Extract data slices only once
    data_shear = np.array(dataset.variables['qshear'])[:pre_ri_time, bottom:, 0, :]
    data_buoy = np.array(dataset.variables['qbuoy'])[:pre_ri_time, bottom:, 0, :]
    data_diss = np.array(dataset.variables['qdiss'])[:pre_ri_time, bottom:, 0, :]
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
    data_shear = np.array(dataset.variables['qshear'])[pre_ri_time:post_ri_time, bottom:, 0, :]
    data_buoy = np.array(dataset.variables['qbuoy'])[pre_ri_time:post_ri_time, bottom:, 0, :]
    data_diss = np.array(dataset.variables['qdiss'])[pre_ri_time:post_ri_time, bottom:, 0, :]
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
    data_shear = np.array(dataset.variables['qshear'])[post_ri_time:, bottom:, 0, :]
    data_buoy = np.array(dataset.variables['qbuoy'])[post_ri_time:, bottom:, 0, :]
    data_diss = np.array(dataset.variables['qdiss'])[post_ri_time:, bottom:, 0, :]
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
