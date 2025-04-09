# coding=utf-8
"""
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm


# Function to plot data
def plot_and_save(data_func: np.array, title: str, filename: str, xticks: np.array,
                  cmap: str, cbar_label='TKE per Second ($m^2/s^3$)') -> None:
    """
    Function for plotting arrays.
    :param xticks: X-axis ticks
    :param data_func: Array to be plotted.
    :param title: Title of plot.
    :param filename: Filename of plot being saved.
    :param cmap: Colormap of plot.
    :param cbar_label: Label of the colorbar.
    :return: Nothing.
    """
    plt.contourf(data_func, cmap=cmap, extent=(np.min(x), np.max(x), np.min(z), np.max(z)))
    plt.xlabel(r'Distance from TC Center (km)')
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(xticks)
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
    output_dir = "//uahdata/rstor/aes655_project/sep_by_intensity_phase/"

    # Load dataset
    with netCDF4.Dataset(path) as ds:
        z = np.asarray(ds.variables['lev'])
        x = np.asarray(ds.variables['lon'])
        mtime = (np.asarray(ds.variables['mtime']) / 3600)[:, 0, 0]
        pre_ri_time = np.where(mtime == 10)[0][0]
        post_ri_time = np.where(mtime == 71)[0][0]
        bottom = np.searchsorted(z, 10) - 1  # More efficient than np.where(z > 10)[0][0]
        top = np.searchsorted(z, 20) + 1
        z = z[bottom: top]

        # Extract and slice data
        variables = ['qshear', 'qbuoy', 'qdiss', 'dqke', 'qke_adv', 'qke', 'qwt']
        preri_data = {var: np.asarray(ds.variables[var])[:pre_ri_time, bottom: top, 0, :] for var in variables}
        ri_data = {var: np.asarray(ds.variables[var])[pre_ri_time:post_ri_time, bottom: top, 0, :] for var in variables}
        postri_data = {var: np.asarray(ds.variables[var])[post_ri_time:, bottom: top, 0, :] for var in variables}

    # Compute average values
    preri_avg_data = {key: np.nanmean(val, axis=0) for key, val in preri_data.items()}
    ri_avg_data = {key: np.nanmean(val, axis=0) for key, val in ri_data.items()}
    postri_avg_data = {key: np.nanmean(val, axis=0) for key, val in postri_data.items()}

    # Plot averaged data
    plot_and_save(preri_avg_data['qshear'], 'Average TKE Shear Production Before RI', 'preri_shear_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(preri_avg_data['qbuoy'], 'Average TKE Buoyancy Production Before RI', 'preri_buoy_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(preri_avg_data['qdiss'], 'Average TKE Dissipation Before RI', 'preri_diss_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(preri_avg_data['dqke'], 'Average TKE Change Before RI', 'preri_change_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(preri_avg_data['qke_adv'], 'Average TKE Advection Before RI', 'preri_adv_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(preri_avg_data['qwt'], 'Average TKE Vertical Transport Before RI', 'preri_vert_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(preri_avg_data['qke'], 'Average TKE Before RI', 'preri_ke_avg.png', cbar_label='TKE ($m^2/s^2$)',
                  cmap='turbo', xticks=np.arange(0, 301, 50))

    # Plot averaged data
    plot_and_save(ri_avg_data['qshear'], 'Average TKE Shear Production During RI', 'ri_shear_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(ri_avg_data['qbuoy'], 'Average TKE Buoyancy Production During RI', 'ri_buoy_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(ri_avg_data['qdiss'], 'Average TKE Dissipation During RI', 'ri_diss_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(ri_avg_data['dqke'], 'Average TKE Change During RI', 'ri_change_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(ri_avg_data['qke_adv'], 'Average TKE Advection During RI', 'ri_adv_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(ri_avg_data['qwt'], 'Average TKE Vertical Transport During RI', 'ri_vert_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(ri_avg_data['qke'], 'Average TKE During RI', 'ri_ke_avg.png', cbar_label='TKE ($m^2/s^2$)',
                  cmap='turbo', xticks=np.arange(0, 301, 50))

    # Plot averaged data
    plot_and_save(postri_avg_data['qshear'], 'Average TKE Shear Production After RI', 'postri_shear_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(postri_avg_data['qbuoy'], 'Average TKE Buoyancy Production After RI', 'postri_buoy_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(postri_avg_data['qdiss'], 'Average TKE Dissipation After RI', 'postri_diss_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(postri_avg_data['dqke'], 'Average TKE Change After RI', 'postri_change_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(postri_avg_data['qke_adv'], 'Average TKE Advection After RI', 'postri_adv_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(postri_avg_data['qwt'], 'Average TKE Vertical Transport After RI', 'postri_vert_avg.png',
                  cmap='turbo', xticks=np.arange(0, 301, 50))
    plot_and_save(postri_avg_data['qke'], 'Average TKE After RI', 'postri_ke_avg.png', cbar_label='TKE ($m^2/s^2$)',
                  cmap='turbo', xticks=np.arange(0, 301, 50))

    # Create colormap for dominant TKE production term
    colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF']  # Red, Green, Blue, Yellow, Magenta
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

    preri_dominant_term = compute_dominant_term(preri_avg_data['qshear'], preri_avg_data['qbuoy'],
                                                preri_avg_data['qke_adv'], preri_avg_data['qdiss'],
                                                preri_avg_data['qwt'])
    ri_dominant_term = compute_dominant_term(ri_avg_data['qshear'], ri_avg_data['qbuoy'],
                                             ri_avg_data['qke_adv'], ri_avg_data['qdiss'],
                                             preri_avg_data['qwt'])
    postri_dominant_term = compute_dominant_term(postri_avg_data['qshear'], postri_avg_data['qbuoy'],
                                                 postri_avg_data['qke_adv'], postri_avg_data['qdiss'],
                                                 preri_avg_data['qwt'])

    # Plot dominant TKE term
    plt.figure(figsize=(10, 6))
    plt.imshow(preri_dominant_term, cmap=custom_cmap, aspect='auto',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
               vmin=1, vmax=5)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Distance from TC Center (km)')
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.title('Dominant TKE Production Term Before RI')
    cbar = plt.colorbar()
    cbar.set_ticks([1, 2, 3, 4, 5])
    cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.savefig(os.path.join(output_dir, 'preri_max_occurrence.png'))
    plt.close()
    # Plot dominant TKE term
    plt.imshow(ri_dominant_term, cmap=custom_cmap, aspect='auto', extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
               vmin=1, vmax=5)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Distance from TC Center (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.ylabel('Height (km)')
    plt.title('Dominant TKE Production Term During RI')
    cbar = plt.colorbar()
    cbar.set_ticks([1, 2, 3, 4, 5])
    cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.savefig(os.path.join(output_dir, 'ri_max_occurrence.png'))
    plt.close()
    # Plot dominant TKE term
    plt.imshow(postri_dominant_term, cmap=custom_cmap, aspect='auto',
               extent=(np.min(x), np.max(x), np.max(z), np.min(z)),
               vmin=1, vmax=5)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Distance from TC Center (km)')
    plt.ylabel('Height (km)')
    plt.yticks(ticks=np.arange(10, 21, 2))
    plt.xticks(ticks=np.arange(0, 301, 50))
    plt.title('Dominant TKE Production Term After RI')
    cbar = plt.colorbar()
    cbar.set_ticks([1, 2, 3, 4, 5])
    cbar.set_ticklabels(['Shear', 'Buoyancy', 'Advection', 'Dissipation', 'Vert. Transport'])
    plt.savefig(os.path.join(output_dir, 'postri_max_occurrence.png'))
    plt.close()
