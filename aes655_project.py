# coding=utf-8
"""
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


# Function to plot data
def plot_and_save(data_func: np.array, title: str, filename: str, cmap='rainbow',
                  cbar_label='TKE per Second ($m^2/s^3$)') -> None:
    """
    Function for plotting arrays.
    :param data_func: Array to be plotted.
    :param title: Title of plot.
    :param filename: Filename of plot being saved.
    :param cmap: Colormap of plot.
    :param cbar_label: Label of the colorbar.
    :return: Nothing.
    """
    plt.imshow(data_func, aspect='auto', cmap=cmap, extent=(np.min(x), np.max(x), np.max(z), np.min(z)))
    plt.gca().invert_yaxis()
    plt.xlabel('Range (km)')
    plt.ylabel('Height (km)')
    plt.title(title)
    plt.colorbar(label=cbar_label)
    plt.savefig(os.path.join(output_dir, filename))
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


if __name__ == "__main__":
    # Define the file path
    path = "//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc"
    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Average_Data/TKE/"

    # Load dataset
    with netCDF4.Dataset(path) as ds:
        z = np.array(ds.variables['lev'])
        x = np.array(ds.variables['lon'])
        bottom = np.where(z > 10)[0][0]
        z = z[bottom:]

        # Extract and slice data
        variables = ['qshear', 'qbuoy', 'qdiss', 'dqke', 'qke_adv', 'qke', 'qwt']
        data = {var: np.array(ds.variables[var])[:, bottom:, 0, :] for var in variables}

    # Compute average values
    avg_data = {key: np.nanmean(val, axis=0) for key, val in data.items()}

    # Plot averaged data
    plot_and_save(avg_data['qshear'], 'Average TKE Shear Production', 'shear_avg.png')
    plot_and_save(avg_data['qbuoy'], 'Average TKE Buoyancy Production', 'buoy_avg.png')
    plot_and_save(avg_data['qdiss'], 'Average TKE Dissipation', 'diss_avg.png')
    plot_and_save(avg_data['dqke'], 'Average TKE Change', 'change_avg.png')
    plot_and_save(avg_data['qke_adv'], 'Average TKE Advection', 'adv_avg.png')
    plot_and_save(avg_data['qwt'], 'Average TKE Vertical Transport', 'vert_avg.png')
    plot_and_save(avg_data['qke'], 'Average TKE', 'ke_avg.png', cbar_label='TKE ($m^2/s^2$)')

    # Create colormap for dominant TKE production term
    colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00']  # Red, Green, Blue, Yellow
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

    dominant_term = compute_dominant_term(avg_data['qshear'], avg_data['qbuoy'], avg_data['qke_adv'], avg_data['qdiss'])

    # Plot dominant TKE term
    plt.imshow(dominant_term, aspect='auto', cmap=custom_cmap, extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
               vmin=1, vmax=4)
    plt.gca().invert_yaxis()
    plt.xlabel('Range (km)')
    plt.ylabel('Height (km)')
    plt.title('Dominant TKE Production Term')
    cbar = plt.colorbar()
    cbar.set_ticks([1, 2, 3, 4])
    cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation'])
    plt.savefig(os.path.join(output_dir, 'max_occurrence.png'))
    plt.close()

    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/"
    # Loop through timesteps
    hourly_dirs = {var: os.path.join(output_dir, f'{var}_hourly') for var in variables}
    for directory in hourly_dirs.values():
        os.makedirs(directory, exist_ok=True)

    for i in range(data['qshear'].shape[0]):
        for var in variables:
            plot_and_save(data[var][i], f'TKE {var} at {i} Timestep', os.path.join(hourly_dirs[var], f'{var}_{i}.png'))

        dominant_term = compute_dominant_term(data['qshear'][i], data['qbuoy'][i], data['qke_adv'][i], data['qdiss'][i])
        plt.imshow(dominant_term, aspect='auto', cmap=custom_cmap, extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
                   vmin=1, vmax=4)
        plt.gca().invert_yaxis()
        plt.xlabel('Range (km)')
        plt.ylabel('Height (km)')
        plt.title(f'Dominant TKE Production Term at {i} Timestep')
        cbar = plt.colorbar()
        cbar.set_ticks([1, 2, 3, 4])
        cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation'])
        plt.savefig(os.path.join(output_dir, 'occ_hourly', f'max_occurrence_{i}.png'))
        plt.close()
