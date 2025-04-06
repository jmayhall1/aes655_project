# coding=utf-8
"""
@author: John Mark Mayhall
"""
import collections
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
    plt.yticks(ticks=np.arange(10, 21, 2))
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
        bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
        top = np.searchsorted(z, 20) + 1
        z = z[bottom:top]

        # Extract and slice data
        variables = ['qshear', 'qbuoy', 'qdiss', 'dqke', 'qke_adv', 'qke', 'qwt']
        data = {var: np.array(ds.variables[var])[:, bottom:top, 0, :] for var in variables}

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
    colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF']  # Red, Green, Blue, Yellow, Magenta
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

    dominant_term = compute_dominant_term(avg_data['qshear'], avg_data['qbuoy'], avg_data['qke_adv'],
                                          avg_data['qdiss'], avg_data['qwt'])

    # Plot dominant TKE term
    plt.imshow(dominant_term, aspect='auto', cmap=custom_cmap, extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
               vmin=1, vmax=5)
    plt.gca().invert_yaxis()
    plt.xlabel('Range (km)')
    plt.ylabel('Height (km)')
    plt.title('Dominant TKE Production Term')
    plt.yticks(ticks=np.arange(10, 21, 2))
    cbar = plt.colorbar()
    cbar.set_ticks([1, 2, 3, 4, 5])
    cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.savefig(os.path.join(output_dir, 'max_occurrence.png'))
    plt.close()

    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/"
    # Loop through timesteps
    hourly_dirs = {var: os.path.join(output_dir, f'{var}_hourly') for var in variables}
    for directory in hourly_dirs.values():
        os.makedirs(directory, exist_ok=True)

    dom_term_list = []
    shear_list, buoy_list, adv_list, diss_list, vert_list = [], [], [], [], []
    for i in range(data['qshear'].shape[0]):
        for var in variables:
            plot_and_save(data[var][i], f'TKE {var} at {i} Timestep', os.path.join(hourly_dirs[var], f'{var}_{i}.png'))

        dominant_term = compute_dominant_term(data['qshear'][i], data['qbuoy'][i], data['qke_adv'][i],
                                              data['qdiss'][i], data['qwt'][i])
        dom_term_list.append(np.argmax(np.bincount(dominant_term.flatten())))
        values = collections.Counter(dominant_term.flatten())
        shear_list.append(values[1])
        buoy_list.append(values[2])
        adv_list.append(values[3])
        diss_list.append(values[4])
        vert_list.append(values[5])
        plt.imshow(dominant_term, aspect='auto', cmap=custom_cmap, extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
                   vmin=1, vmax=5)
        plt.gca().invert_yaxis()
        plt.xlabel('Range (km)')
        plt.yticks(ticks=np.arange(10, 21, 2))
        plt.ylabel('Height (km)')
        plt.title(f'Dominant TKE Production Term at {i} Timestep')
        cbar = plt.colorbar()
        cbar.set_ticks([1, 2, 3, 4, 5])
        cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
        plt.savefig(os.path.join(output_dir, 'occ_hourly', f'max_occurrence_{i}.png'))
        plt.close()

    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Line_Data/"
    time = np.arange(data['qshear'].shape[0])
    # Plot max term over time
    plt.figure(figsize=(10, 6))
    plt.plot(time, dom_term_list, marker='o', linestyle='-')
    plt.ylabel('Dominant Term')
    plt.xlabel('Time Step')
    plt.yticks(ticks=np.arange(1, 6), labels=['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.title('Dominant Term vs Time')
    plt.savefig(os.path.join(output_dir, 'max_occ_over_time.png'))
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(time, shear_list, marker='o', linestyle='-', label='Shear')
    plt.plot(time, buoy_list, marker='o', linestyle='-', label='Buoyancy')
    plt.plot(time, adv_list, marker='o', linestyle='-', label='Advection')
    plt.plot(time, diss_list, marker='o', linestyle='-', label='Dissapation')
    plt.plot(time, vert_list, marker='o', linestyle='-', label='Vert. Transport')
    plt.ylabel('Pixel Count')
    plt.xlabel('Time Step')
    plt.title('Pixel Count for Terms vs Time')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'term_count.png'))
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(time, shear_list, marker='o', linestyle='-', label='Shear')
    plt.plot(time, buoy_list, marker='o', linestyle='-', label='Buoyancy')
    plt.plot(time, adv_list, marker='o', linestyle='-', label='Advection')
    plt.plot(time, diss_list, marker='o', linestyle='-', label='Dissapation')
    plt.plot(time, vert_list, marker='o', linestyle='-', label='Vert. Transport')
    plt.ylabel('Pixel Count')
    plt.ylim((0, 500))
    plt.xlabel('Time Step')
    plt.title('Pixel Count for Terms vs Time')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'term_count_zoomed.png'))
    plt.close()
