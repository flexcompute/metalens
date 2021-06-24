import numpy as np
import matplotlib.pylab as plt
import matplotlib

import tidy3d as td

# 1 nanometer in units of microns (for conversion)
nm = 1e-3

# free space central wavelength
wavelength = 660 * nm

# desired numerical aperture
NA = 0.8

# shape parameters of metalens unit cell (um) (refer to image above and see paper for details)
rect_width = 85 * nm
rect_length = 410 * nm
lens_thick = 600 * nm
cell_length = 430 * nm

# space between PML and bottom of metalens, also space between top of metalens and PML
spacing = 1 * wavelength

# side length of entire metalens (um)
side_length = 21 / 1.5 * wavelength

# Number of unit cells in each x and y direction (NxN grid)
N = int(side_length / cell_length)
print(f'for diameter of {(side_length / wavelength):.1f} wavelengths, have {N} cells per side')
print(f'full metalens has area of {(side_length / wavelength)**2:.1f} wavelengths^2 and {N*N} total cells')

# Define material properties at 600 nm
n_TiO2 = 2.40
n_SiO2 = 1.46
air = td.Medium(epsilon=1.0)
SiO2 = td.Medium(epsilon=n_SiO2**2)
TiO2 = td.Medium(epsilon=n_TiO2**2)

# resolution control
grids_per_wavelength = 18

# Number of PML layers to use along z direction
npml = 15


# grid size (um)
dl = wavelength / grids_per_wavelength

# using the wavelength in microns, one can use td.C_0 (um/s) to get frequency in Hz
freq_0 = td.C_0 / wavelength
df = freq_0 / 40

# Define PML layers, for this we have no PML in x, y but `npml` cells in z
pml_layers = [npml, npml, npml]

# Compute the domain size in x, y (note: round down from side_length)
length_xy = N * cell_length

# focal length given diameter and numerical aperture
focal_length = length_xy / 2 / NA * np.sqrt(1 - NA**2)
print(f'with NA of {NA}, focal length is {(focal_length / wavelength):.2f} wavelengths or {focal_length:.2f} um')

# Function describing the theoretical best angle of each box at position (x,y).  see paper for details
def theta(x, y):
    return np.pi / wavelength * (focal_length - np.sqrt(x**2 + y**2 + focal_length**2))

length_z = spacing + lens_thick + 1.1 * focal_length + spacing

# construct simulation size array
sim_size = np.array([length_xy, length_xy, length_z])

# define coordinates of each unit cell
centers_x = cell_length * np.arange(N) - length_xy/2. + cell_length/2.
centers_y = cell_length * np.arange(N) - length_xy/2. + cell_length/2.
center_z = -length_z/2. + spacing + lens_thick/2.

# convenience function to make an angled box at each x,y location using polyslab.
# For more details see, https://simulation.cloud/docs/html/generated/tidy3dclient.PolySlab.html
def angled_box(x, y, angle):
    """ make a box of size (L, W, H) centered at `(x, y)` at `angle` from x axis"""

    # x, y vertices of box of size (rect_length, rect_width) centered at the origin
    vertices_origin = np.array([[+rect_length/2, +rect_width/2],
                                [-rect_length/2, +rect_width/2],
                                [-rect_length/2, -rect_width/2],
                                [+rect_length/2, -rect_width/2]])
    
    # 2x2 rotation matrix angle `angle` with respect to x axis
    rotation_matrix = np.array([[+np.cos(angle), -np.sin(angle)],
                                [+np.sin(angle), +np.cos(angle)]])

    # rotate the origin vertices by this angle
    vertices_rotated = vertices_origin @ rotation_matrix
    
    # shift the rotated vertices to be centered at (x, y)
    vertices = vertices_rotated + np.array([x, y])

    # create a tidy3D PolySlab with these rotated and shifted vertices and thickness `H`
    return td.PolySlab(
            vertices,
            z_cent=center_z,
            z_size=lens_thick,
            material=TiO2,
            name=str(np.random.random())
        )

starting_angles = np.zeros((centers_x.size, centers_y.size), dtype=float)

# loop through the coordinates and add all unit cells to geometry list
for i, x in enumerate(centers_x):
    for j, y in enumerate(centers_y):
        angle = theta(x, y)
        starting_angles[i, j] = angle

def create_geometry(angles):
    substrate = td.Box(
       center=[0, 0, 0.0],
       size=[td.inf, td.inf, 1.0],
       material=None)

    geometry = [substrate]
    for i, x in enumerate(centers_x):
        for j, y in enumerate(centers_y):
            theta_ij = angles[i, j]
            geometry.append(angled_box(x, y, theta_ij))
    return geometry

def plot_angles(angles, facecolor='grey'):
    geometry = create_geometry(angles)
    plot_geometry(geometry, facecolor=facecolor)

def plot_geometry(geometry, facecolor='grey', len_frac=1.0):
    substrate, *metacells = geometry

    def get_patch(tidy_polyslab):
        # see https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.Patch.html#matplotlib.patches.Patch
        xy_vertices = tidy_polyslab.vertices
        return matplotlib.patches.Polygon(
            xy_vertices,
            capstyle='round',
            facecolor=facecolor,
            edgecolor='#3b3736',
            hatch=None,    
            zorder=-np.inf,
        )    

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(1, 1, 1)
    for cell in metacells:
        ax.add_patch(get_patch(cell))
    ax.set_xlim([len_frac*-length_xy/2, len_frac*length_xy/2])
    ax.set_ylim([len_frac*-length_xy/2, len_frac*length_xy/2])
    plt.show()

effs = np.load('effs_N21.npy')
eff_best = np.load('eff_best_N21.npy')
improved_angles = np.load('improved_angles_N21.npy')
delta_angles = improved_angles - starting_angles

# plot_angles(starting_angles, facecolor='red')
# plot_angles(improved_angles, facecolor='green')

# imshow
# plt.imshow(delta_angles / 2 / np.pi, cmap='RdBu')
# plt.colorbar()
# plt.xlabel('x index')
# plt.ylabel('y index')
# plt.title('change in angle after optimization $(2 \pi)$')
# plt.show()

# scatter
ipoints, jpoints = np.meshgrid(np.arange(21), np.arange(21))
ipoints = ipoints.flatten()
jpoints = jpoints.flatten()
deltapoints = delta_angles.flatten()

pos = (deltapoints > 0)
neg = ~pos
scale = 200
plt.scatter(ipoints[pos], jpoints[pos], s=scale*np.abs(deltapoints[pos]), color='blue', edgecolor='grey')
plt.scatter(ipoints[neg], jpoints[neg], s=scale*np.abs(deltapoints[neg]), color='red', edgecolor='grey')
plt.xticks(np.arange(21))
plt.yticks(np.arange(21))
plt.xlabel('x index')
plt.ylabel('y index')
plt.title('change in angle after optimization')
plt.show()


# N = 21
eff_start = effs[0]

# N = 10
# eff_start = 0.517

t_solve = 4.06333  # from log

# calculations
num_sims = 1 + np.arange(len(effs))
effs_perc = 100 * np.array(effs)
effs_improved = 100 * (np.array(effs) - eff_start) / eff_start

# f, ((ax1, ax2),
#     (ax3, ax4)) = plt.subplots(2, 2, tight_layout=True, figsize=(10, 6))

# ax1.plot(num_sims, effs_perc)
# ax1.set_xlabel('number of simulations')
# ax1.set_ylabel('efficiency (%)')

# ax2.plot(t_solve/60 * num_sims, effs_perc)
# ax2.set_xlabel('solver time (min)')
# ax2.set_ylabel('efficiency (%)')

# ax3.plot(effs_improved)
# ax3.set_xlabel('number of simulations')
# ax3.set_ylabel('efficiency improvement (%)')

# ax4.plot(t_solve/60 * num_sims, effs_improved)
# ax4.set_xlabel('solver time (min)')
# ax4.set_ylabel('efficiency improvement (%)')

# plt.show()