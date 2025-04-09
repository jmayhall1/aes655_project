# coding=utf-8
"""
Last Edited: 04/09/2025
@author: John Mark Mayhall
"""
import collections
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm


# Function to plot data
def plot_and_save(data_func: np.array, title: str, filename: str, levels: np.array, xticks: np.array, cticks: np.array,
                  cmap: str, cbar_label='TKE per Second ($m^2/s^3$)', norm=None) -> None:
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
    :return: Nothing.
    """
    plt.figure(figsize=(10, 6))
    plt.contourf(data_func, cmap=cmap, extent=(np.min(x), np.max(x), np.min(z), np.max(z)),
                 levels=levels, norm=norm)
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


if __name__ == "__main__":
    # Define the file path
    path = "//uahdata/rstor/aes655_project/cm1out_azimavg_s.nc"
    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Average_Data/TKE/"

    # Load dataset
    with netCDF4.Dataset(path) as ds:
        z = np.asarray(ds.variables['lev'])
        x = np.asarray(ds.variables['lon'])
        bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
        top = np.searchsorted(z, 20) + 1
        z = z[bottom:top]

        # Extract and slice data
        variables = ['qshear', 'qbuoy', 'qdiss', 'dqke', 'qke_adv', 'qke', 'qwt']
        data = {var: np.asarray(ds.variables[var])[:, bottom:top, 0, :] for var in variables}

    # Compute average values
    avg_data = {key: np.nanmean(val, axis=0) for key, val in data.items()}

    # Plot averaged data
    plot_and_save(avg_data['qshear'], 'Average TKE Shear Production', 'shear_avg.png', cmap='turbo',
                  levels=np.linspace(0.0, 0.04, 21), cticks=np.arange(0, 0.04, 0.005),
                  xticks=np.arange(0, 301, 50))
    plot_and_save(avg_data['qbuoy'], 'Average TKE Buoyancy Production', 'buoy_avg.png', cmap='turbo',
                  levels=np.linspace(-0.005, 0.005, 22), cticks=np.arange(-0.005, 0.005, 0.001),
                  xticks=np.arange(0, 301, 50), norm=TwoSlopeNorm(vmin=-0.005, vcenter=0, vmax=0.005))
    plot_and_save(avg_data['qdiss'], 'Average TKE Dissipation', 'diss_avg.png', cmap='turbo',
                  levels=np.linspace(-0.025, 0.0001, 24), cticks=np.arange(-0.025, 0.01, 0.005),
                  xticks=np.arange(0, 301, 50))
    plot_and_save(avg_data['dqke'], 'Average TKE Change', 'change_avg.png', cmap='turbo',
                  levels=np.linspace(-0.001, 0.0031, 19), cticks=np.arange(-0.001, 0.0031, 0.0005),
                  xticks=np.arange(0, 301, 50))
    plot_and_save(avg_data['qke_adv'], 'Average TKE Advection', 'adv_avg.png', cmap='turbo',
                  levels=np.linspace(-0.01, 0.0051, 19), cticks=np.arange(-0.01, 0.0051, 0.005),
                  xticks=np.arange(0, 301, 50), norm=TwoSlopeNorm(vmin=-0.01, vcenter=0, vmax=0.005))
    plot_and_save(avg_data['qwt'], 'Average TKE Vertical Transport', 'vert_avg.png', cmap='turbo',
                  levels=np.linspace(-0.01, 0.005, 24), cticks=np.arange(-0.01, 0.0051, 0.001),
                  xticks=np.arange(0, 301, 50), norm=TwoSlopeNorm(vmin=-0.01, vcenter=0, vmax=0.005))
    plot_and_save(avg_data['qke'], 'Average TKE', 'ke_avg.png', cbar_label='TKE ($m^2/s^2$)', cmap='turbo',
                  levels=np.linspace(-0.001, 10, 22), cticks=np.arange(0, 10, 2),
                  xticks=np.arange(0, 301, 50))

    # Create colormap for dominant TKE production term
    colors_lst = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF']  # Red, Green, Blue, Yellow, Magenta
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors_lst)

    dominant_term = compute_dominant_term(avg_data['qshear'], avg_data['qbuoy'], avg_data['qke_adv'],
                                          avg_data['qdiss'], avg_data['qwt'])

    # Plot dominant TKE term
    plt.figure(figsize=(10, 6))
    plt.imshow(dominant_term, aspect='auto', cmap=custom_cmap, extent=(np.min(x), np.max(x), np.max(z),
                                                                       np.min(z)), vmin=1, vmax=5)
    plt.xlabel(r'Distance from TC Center (km)')
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.title('Dominant TKE Production Term')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    cbar = plt.colorbar()
    cbar.set_ticks([1, 2, 3, 4, 5])
    cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.savefig(os.path.join(output_dir, 'max_occurrence.png'), dpi=300)
    plt.close()

    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Hourly_Data/"
    # Loop through timesteps
    hourly_dirs = {var: os.path.join(output_dir, f'{var}_hourly') for var in variables}
    for directory in hourly_dirs.values():
        os.makedirs(directory, exist_ok=True)

    dom_term_list = []
    shear_list, buoy_list, adv_list, diss_list, vert_list = [], [], [], [], []
    for i in range(data['qshear'].shape[0]):
        plot_and_save(data['qshear'][i], f'TKE Shear Production at Hour {i}',
                      os.path.join(hourly_dirs['qshear'], f'qshear_{i}.png'), cmap='turbo',
                      levels=np.linspace(0.001, 0.2, 21), cticks=np.arange(0, 0.21, 0.02),
                      xticks=np.arange(0, 301, 50))
        plot_and_save(data['qbuoy'][i], f'TKE Buoyancy Production at Hour {i}',
                      os.path.join(hourly_dirs['qbuoy'], f'qbuoy_{i}.png'), cmap='turbo',
                      levels=np.linspace(-0.015, 0.0151, 22), cticks=np.arange(-0.015, 0.0151, 0.005),
                      xticks=np.arange(0, 301, 50), norm=TwoSlopeNorm(vmin=-0.02, vcenter=0, vmax=0.02))
        plot_and_save(data['qdiss'][i], f'TKE Dissipation at Hour {i}',
                      os.path.join(hourly_dirs['qdiss'], f'qdiss_{i}.png'), cmap='turbo',
                      levels=np.linspace(-0.12, -0.0001, 24), cticks=np.arange(-0.12, 0.01, 0.01),
                      xticks=np.arange(0, 301, 50))
        plot_and_save(data['dqke'][i], f'TKE Change at Hour {i}',
                      os.path.join(hourly_dirs['dqke'], f'dqke_{i}.png'), cmap='turbo',
                      levels=np.linspace(-0.01, 0.021, 19), cticks=np.arange(-0.01, 0.021, 0.005),
                      xticks=np.arange(0, 301, 50))
        plot_and_save(data['qke_adv'][i], f'TKE Advection at Hour {i}',
                      os.path.join(hourly_dirs['qke_adv'], f'qke_adv_{i}.png'), cmap='turbo',
                      levels=np.linspace(-0.06, 0.02, 19), cticks=np.arange(-0.06, 0.021, 0.01),
                      xticks=np.arange(0, 301, 50), norm=TwoSlopeNorm(vmin=-0.06, vcenter=0, vmax=0.02))
        plot_and_save(data['qwt'][i], f'TKE Vertical Transport at Hour {i}',
                      os.path.join(hourly_dirs['qwt'], f'qwt_{i}.png'), cmap='turbo',
                      levels=np.linspace(-0.06, 0.05, 24), cticks=np.arange(-0.06, 0.051, 0.01),
                      xticks=np.arange(0, 301, 50), norm=TwoSlopeNorm(vmin=-0.06, vcenter=0, vmax=0.05))
        plot_and_save(data['qke'][i], f'TKE at Hour {i}',
                      os.path.join(hourly_dirs['qke'], f'qke_{i}.png'), cbar_label='TKE ($m^2/s^2$)', cmap='turbo',
                      levels=np.linspace(0.001, 44, 22), cticks=np.arange(0, 46, 5),
                      xticks=np.arange(0, 301, 50))

        dominant_term = compute_dominant_term(data['qshear'][i], data['qbuoy'][i], data['qke_adv'][i],
                                              data['qdiss'][i], data['qwt'][i])
        dom_term_list.append(np.argmax(np.bincount(dominant_term.flatten())))
        values = collections.Counter(dominant_term.flatten())
        shear_list.append(values[1])
        buoy_list.append(values[2])
        adv_list.append(values[3])
        diss_list.append(values[4])
        vert_list.append(values[5])
        plt.figure(figsize=(10, 6))
        plt.imshow(dominant_term, aspect='auto', cmap=custom_cmap, extent=(np.min(x), np.max(x), np.max(z),
                                                                           np.min(z)), vmin=1, vmax=5)
        plt.xlabel(r'Distance from TC Center (km)')
        plt.gca().invert_yaxis()
        plt.yticks(ticks=np.arange(10, 21, 2))
        plt.xticks(np.arange(0, 301, 50))
        plt.ylabel('Height (km)')
        plt.title(f'Dominant TKE Production Term at Hour {i}')
        cbar = plt.colorbar()
        cbar.set_ticks([1, 2, 3, 4, 5])
        cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
        plt.savefig(os.path.join(output_dir, 'occ_hourly', f'max_occurrence_{i}.png'), dpi=300)
        plt.close()

    output_dir = "//uahdata/rstor/aes655_project/not_sep_by_intensity_phase/Line_Data/"
    time = np.arange(data['qshear'].shape[0])
    # Plot max term over time
    plt.figure(figsize=(10, 6))
    plt.plot(time, dom_term_list, marker='o', linestyle='-')
    plt.ylabel('Dominant Term')
    plt.xlabel('Hour')
    plt.yticks(ticks=np.arange(1, 6), labels=['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.title('Dominant Term vs Time')
    plt.savefig(os.path.join(output_dir, 'max_occ_over_time.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(time, shear_list, marker='o', linestyle='-', label='Shear')
    plt.plot(time, buoy_list, marker='o', linestyle='-', label='Buoyancy')
    plt.plot(time, adv_list, marker='o', linestyle='-', label='Advection')
    plt.plot(time, diss_list, marker='o', linestyle='-', label='Dissipation')
    plt.plot(time, vert_list, marker='o', linestyle='-', label='Vert. Transport')
    plt.ylabel('# of Pixels that are Dominated by a TKE Term')
    plt.xlabel('Hour')
    plt.title('Pixel Count for Terms vs Time')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'term_count.png'), dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(time, shear_list, marker='o', linestyle='-', label='Shear')
    plt.plot(time, buoy_list, marker='o', linestyle='-', label='Buoyancy')
    plt.plot(time, adv_list, marker='o', linestyle='-', label='Advection')
    plt.plot(time, diss_list, marker='o', linestyle='-', label='Dissipation')
    plt.plot(time, vert_list, marker='o', linestyle='-', label='Vert. Transport')
    plt.ylabel('# of Pixels that are Dominated by a TKE Term')
    plt.ylim((0, 500))
    plt.xlabel('Hour')
    plt.title('Pixel Count for Terms vs Time Excluding Higher Counts')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'term_count_zoomed.png'), dpi=300)
    plt.close()
