# coding=utf-8
"""
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np

# The ile path
path = '//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc'

# Load netCDF dataset once
dataset = netCDF4.Dataset(path)

# Read variables
z = np.array(dataset.variables['lev'])
x = np.array(dataset.variables['lon'])
bottom = np.searchsorted(z, 10)  # More efficient than np.where(z > 10)[0][0]
z = z[bottom:]

# Extract data slices only once
data_shear = np.array(dataset.variables['qshear'])[:, bottom:, 0, :]
data_buoy = np.array(dataset.variables['qbuoy'])[:, bottom:, 0, :]
data_diss = np.array(dataset.variables['qdiss'])[:, bottom:, 0, :]


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


# Compute normalized values
buoy_shear = safe_divide(data_buoy, data_shear)
diss_shear = safe_divide(data_diss, data_shear)
buoydiss_shear = safe_divide(data_buoy + data_diss, data_shear)

# Compute averages
buoy_avg_shear = np.nanmean(buoy_shear, axis=0)
diss_avg_shear = np.nanmean(diss_shear, axis=0)
buoydiss_avg_shear = np.nanmean(buoydiss_shear, axis=0)


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
    plt.savefig(f'//uahdata/rstor/aes655_project/{filename}.png')
    plt.close()


# Save average plots
plot_and_save(buoy_avg_shear, 'Average Buoyancy Production %', 'buoy_avg_percent')
plot_and_save(diss_avg_shear, 'Average Dissipation %', 'diss_avg_percent')
plot_and_save(buoydiss_avg_shear, 'Average Buoyancy Production and Dissipation %', 'buoydiss_avg_percent')

# Ensure output directories exist
os.makedirs('//uahdata/rstor/aes655_project/buoy_percent_hourly', exist_ok=True)
os.makedirs('//uahdata/rstor/aes655_project/diss_percent_hourly', exist_ok=True)
os.makedirs('//uahdata/rstor/aes655_project/buoydiss_percent_hourly', exist_ok=True)

# Generate and save time-step plots
for i in range(data_shear.shape[0]):
    plot_and_save(buoy_shear[i], f'Buoyancy Production % at {i} Timestep',
                  f'buoy_percent_hourly/buoy_percent_{i}')
    plot_and_save(diss_shear[i], f'Dissipation % at {i} Timestep',
                  f'diss_percent_hourly/diss_percent_{i}')
    plot_and_save(buoydiss_shear[i], f'Buoyancy Production and Dissipation % at {i} Timestep',
                  f'buoydiss_percent_hourly/buoydiss_percent_{i}')
