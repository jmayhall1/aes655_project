# coding=utf-8
"""
@author: John Mark Mayhall
"""
import os

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import netCDF4
import numpy as np

# Define file paths
base_path = '//uahdata/rstor/aes655_project/'
output_path = os.path.join(base_path, 'Rb_hourly')
os.makedirs(output_path, exist_ok=True)

# Load NetCDF data efficiently
ncfile = netCDF4.Dataset(os.path.join(base_path, 'cm1out_azimavg_s.nc'))
z = np.array(ncfile.variables['lev'])
x = np.array(ncfile.variables['lon'])
u = np.array(ncfile.variables['u'])[:, :, 0, :]
v = np.array(ncfile.variables['v'])[:, :, 0, :]
thv = np.array(ncfile.variables['thv'])[:, :, 0, :]
ncfile.close()  # Close file after loading data

# Compute vertical differences
delta_z = np.diff(z, prepend=z[0])

# Compute vertical derivatives using MetPy
delta_u = mpcalc.first_derivative(u, axis=1, x=z)
delta_v = mpcalc.first_derivative(v, axis=1, x=z)
delta_thv = mpcalc.first_derivative(thv, axis=1, x=z)

# Find the bottom index where z > 10
bottom = np.searchsorted(z, 10)

# Trim arrays to remove unnecessary levels
z = np.tile(z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
z = np.repeat(z, u.shape[0], axis=0)[:, bottom:, :] * 1000

delta_z = np.tile(delta_z[:, np.newaxis], (1, x.shape[0]))[np.newaxis, :, :]
delta_z = np.repeat(delta_z, u.shape[0], axis=0)[:, bottom:, :] * 1000

u, v, thv = u[:, bottom:, :], v[:, bottom:, :], thv[:, bottom:, :]
delta_u, delta_v, delta_thv = delta_u[:, bottom:, :], delta_v[:, bottom:, :], delta_thv[:, bottom:, :]

# Compute Richardson number
R = (9.81 * delta_z * delta_thv) / (thv * (delta_u ** 2 + delta_v ** 2))
R[np.isinf(R)] = np.nan

# Compute the average Richardson number
R_avg = np.nanmean(R, axis=0)

# Save the average Richardson number plot
plt.figure(figsize=(8, 6))
plt.imshow(R_avg, aspect='auto', cmap='rainbow', vmin=-1, vmax=10,
           extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
plt.gca().invert_yaxis()
plt.ylabel('Height (km)')
plt.xlabel('Range (km)')
plt.title('Average Richardson Number')
plt.colorbar(label=r'$R_b$ (unitless)')
plt.savefig(os.path.join(base_path, 'Rb_avg.png'))
plt.close()

# Initialize tracking arrays
num_timesteps, zdim, xdim = R.shape
minimum_loc = np.zeros((zdim, xdim), dtype=int)
loc025 = np.zeros((zdim, xdim), dtype=int)
loc1 = np.zeros((zdim, xdim), dtype=int)

time = np.arange(num_timesteps)
minimum_values = np.nanmin(R, axis=(1, 2))

# Loop through timesteps for visualization and tracking
for i in range(num_timesteps):
    Ri = R[i, :, :]

    # Identify critical regions
    minimum_loc += (Ri == np.nanmin(Ri)).astype(int)
    loc025 += (Ri < 0.25).astype(int)
    loc1 += (Ri < 1).astype(int)

    # Save individual timestep plot
    plt.figure(figsize=(8, 6))
    plt.imshow(Ri, aspect='auto', cmap='rainbow', vmin=-1, vmax=1,
               extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel('Range (km)')
    plt.title(f'$R_b$ at Timestep {i}')
    plt.colorbar(label=r'$R_b$ (unitless)')
    plt.savefig(os.path.join(output_path, f'Rb_{i}.png'))
    plt.close()

# Plot minimum Richardson number over time
plt.figure(figsize=(8, 6))
plt.plot(time, minimum_values, marker='o', linestyle='-')
plt.ylabel('$R_b$')
plt.xlabel('Time Step')
plt.ylim(0, 2)
plt.title('Minimum $R_b$ vs Time')
plt.savefig(os.path.join(base_path, 'Rb_min.png'))
plt.close()


# Function to save occurrence heatmaps
def save_heatmap(data: np.array, title: str, filename: str, colorbar_label: str) -> None:
    """
    Function to plot occurrence heatmaps.
    :param data: Data to be used for plotting
    :param title: Title of the created plot
    :param filename: Saved plot's filename
    :param colorbar_label: Colorbar Label
    :return: Nothing
    """
    plt.figure(figsize=(8, 6))
    plt.imshow(data, aspect='auto', cmap='rainbow', vmin=0,
               extent=(np.min(x), np.max(x), np.max(z) / 1000, np.min(z) / 1000))
    plt.gca().invert_yaxis()
    plt.ylabel('Height (km)')
    plt.xlabel('Range (km)')
    plt.title(title)
    plt.colorbar(label=colorbar_label)
    plt.savefig(os.path.join(base_path, filename))
    plt.close()


# Save occurrence heatmaps
save_heatmap(minimum_loc, 'Location of Minimum $R_b$', 'minRb_loc.png', r'# of Minimum $R_b$ Occurrences')
save_heatmap(loc025, 'Location of $R_b$ < 0.25', '025Rb_loc.png', r'# of $R_b$ < 0.25 Occurrences')
save_heatmap(loc1, 'Location of $R_b$ < 1', '1Rb_loc.png', r'# of $R_b$ < 1 Occurrences')

