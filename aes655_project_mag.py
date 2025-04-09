# coding=utf-8
"""
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np


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


# Function to generate plots
def plot_and_save(data: np.array, title: str, filename: str) -> None:
    """
    Function for plotting data
    :param data: Array to be plotted
    :param title: Title of plot.
    :param filename: Filename of the plot that is being saved.
    :return: Nothing
    """
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
    plt.savefig(f'//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/{filename}.png')
    plt.close()


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
    z = z[bottom:top]

    # Extract data slices only once
    data_shear = np.asarray(dataset.variables['qshear'])[:, bottom:top, 0, :]
    data_buoy = np.asarray(dataset.variables['qbuoy'])[:, bottom:top, 0, :]
    data_diss = np.asarray(dataset.variables['qdiss'])[:, bottom:top, 0, :]

    # Compute normalized values
    buoy_shear = safe_divide(data_buoy, data_shear)
    buoy_shear[buoy_shear > 110] = 110
    diss_shear = safe_divide(data_diss, data_shear)
    diss_shear[diss_shear > 110] = 110
    buoydiss_shear = safe_divide(data_buoy + data_diss, data_shear)
    buoydiss_shear[buoydiss_shear > 110] = 110

    # Compute averages
    buoy_avg_shear = np.nanmean(buoy_shear, axis=0)
    diss_avg_shear = np.nanmean(diss_shear, axis=0)
    buoydiss_avg_shear = np.nanmean(buoydiss_shear, axis=0)

    # Save average plots
    plot_and_save(buoy_avg_shear, 'Average Buoyancy Production %',
                  'Average_Data/Magnitude/buoy_avg_percent')
    plot_and_save(diss_avg_shear, 'Average Dissipation %', 'Average_Data/Magnitude/diss_avg_percent')
    plot_and_save(buoydiss_avg_shear, 'Average Buoyancy Production and Dissipation %',
                  'Average_Data/Magnitude/buoydiss_avg_percent')

    # Ensure output directories exist
    os.makedirs('//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/buoy_percent_hourly',
                exist_ok=True)
    os.makedirs('//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/diss_percent_hourly',
                exist_ok=True)
    os.makedirs('//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/buoydiss_percent_hourly',
                exist_ok=True)

    # Generate and save time-step plots
    for i in range(data_shear.shape[0]):
        plot_and_save(buoy_shear[i], f'Buoyancy Production % at {i} Hour',
                      f'Hourly_Data/buoy_percent_hourly/buoy_percent_{i}')
        plot_and_save(diss_shear[i], f'Dissipation % at {i} Hour',
                      f'Hourly_Data/diss_percent_hourly/diss_percent_{i}')
        plot_and_save(buoydiss_shear[i], f'Buoyancy Production and Dissipation % at {i} Hour',
                      f'Hourly_Data/buoydiss_percent_hourly/buoydiss_percent_{i}')
