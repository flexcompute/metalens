
# standard python imports
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm

# simulated annealer
from scipy.optimize import minimize

# tidy3D import
import tidy3d as td
import tidy3d.web as web

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
side_length = 10 / 1.5 * wavelength

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


# In[20]:


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
            material=TiO2,
            vertices=vertices,
            z_cent=center_z,
            z_size=lens_thick,
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
       material=SiO2,
       center=[0, 0, -length_z+spacing],
       size=[td.inf, td.inf, length_z])

    geometry = [substrate]
    for i, x in enumerate(centers_x):
        for j, y in enumerate(centers_y):
            theta_ij = angles[i, j]
            geometry.append(angled_box(x, y, theta_ij))
    return geometry

def plot_angles(angles):
    geometry = create_geometry(angles)
    plot_geometry(geometry)

def plot_geometry(geometry, len_frac=1.0):
    substrate, *metacells = geometry

    def get_patch(tidy_polyslab):
        # see https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.Patch.html#matplotlib.patches.Patch
        xy_vertices = tidy_polyslab.vertices
        return matplotlib.patches.Polygon(
            xy_vertices,
            capstyle='round',
            facecolor='grey',
            edgecolor='#3b3736',
            hatch=None,    
            zorder=-np.inf,
        )    

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1)
    for cell in metacells:
        ax.add_patch(get_patch(cell))
    ax.set_xlim([len_frac*-length_xy/2, len_frac*length_xy/2])
    ax.set_ylim([len_frac*-length_xy/2, len_frac*length_xy/2])
    plt.show()

# Bandwidth in Hz
freq_width = freq_0 / 10.

# time dependence of source
gaussian_0 = td.GaussianPulse(freq_0, freq_width, phase=0)
gaussian_i = td.GaussianPulse(freq_0, freq_width, phase=-np.pi/2)

source_x = td.PlaneWave(
            source_time=gaussian_0,
            injection_axis='+z',
            position=-length_z/2 + 2*dl, # just next to PML
            polarization='x')
source_y = td.PlaneWave(
            source_time=gaussian_i,
            injection_axis='+z',
            position=-length_z/2 + 2*dl, # just next to PML
            polarization='y')

# circular polarization Ex + iEy
source_circ = [source_x, source_y]

# Simulation run time past the source decay
run_time = 20. / freq_width

spot_FWHM = 3 * 0.42

monitor_in = td.FreqMonitor(
                center=[0., 0., -length_z/2 + spacing -2*dl],
                size=[length_xy, length_xy, 0],
                freqs=[freq_0],
                store=['E', 'H'],
                name='incident')
monitor_3fwhm = td.FreqMonitor(
                center=[0., 0., -length_z/2 + spacing + lens_thick + focal_length],
                size=[spot_FWHM, spot_FWHM, 0],
                freqs=[freq_0],
                store=['E', 'H'],
                name='focus')
monitors = [monitor_in, monitor_3fwhm]


def run(sim, task_name=''):
    task_id = web.run(sim, task_name=task_name)
    web.monitor_project(task_id)
    web.download_results(task_id)
    sim.load_results('out/monitor_data.hdf5')

def measure_efficiency(sim):
    power_in = sim.flux(monitor_in)[0][0]
    power_focus = sim.flux(monitor_3fwhm)[0][0]
    eff_focus = power_focus / power_in
    return eff_focus

td.logging_level('error')

def compute_efficiency(angles):
    geometry = create_geometry(angles)
    sim = td.Simulation(size=sim_size,
                        mesh_step=[dl, dl, dl],
                        structures=geometry,
                        sources=source_circ,
                        monitors=monitors,
                        run_time=run_time,
                        pml_layers=pml_layers)
    run(sim, task_name='spline_perturb')
    eff_focus = measure_efficiency(sim)
    return eff_focus

# plot_angles(starting_angles)
import time
t = time.time()
eff_design = compute_efficiency(starting_angles)
t_single = time.time() - t
print(f'baseline efficiency = {(eff_design*100):.7f}%')

# amount of perturbation
delta_theta_deg = 5
delta_theta_rad = 2 * np.pi * delta_theta_deg / 360.0

def perturb(i, j, improved_angles, delta_theta):
    perturbed_angles = improved_angles.copy()
    
    perturbed_angles[i, j] += delta_theta
#     perturbed_angles[N//2 + i, j] += delta_theta
#     perturbed_angles[i, N//2 + j] += delta_theta
#     perturbed_angles[N//2 + i, N//2 + j] += delta_theta

    eff_perturbed = compute_efficiency(perturbed_angles)
    return eff_perturbed, perturbed_angles

def line_search(i, j, improved_angles, eff_best, effs):

    print(25*'-')
    print(f'working on i={i+1}/{N}, j={j+1}{N}')

    # do a forward run and determine the sign of the line search
    eff_plus, angles_plus = perturb(i, j, improved_angles, delta_theta=+delta_theta_rad)
    effs.append(eff_best)
    save_progress(improved_angles, eff_best, effs)

    if eff_plus > eff_best:
        # record this change and go forwards
        optimal_step = delta_theta_rad
        eff_best = eff_plus
        improved_angles = angles_plus.copy()
        print(f'+++ efficiency improved to {(eff_best*100):.7f}%: engaging in line search in + direction')
    else:
        optimal_step = -delta_theta_rad

    # line search
    while True:
        eff_perturb, angles_perturb = perturb(i, j, improved_angles, delta_theta=optimal_step)
        effs.append(eff_best)
        save_progress(improved_angles, eff_best, effs)
        if eff_perturb > eff_best:
            print(f'+++ efficiency improved to {(eff_best*100):.7f}%: continuing with line search')
            eff_best = eff_perturb
            improved_angles = angles_perturb.copy()
        else:
            break

    return improved_angles, eff_best
   
def save_progress(improved_angles, eff_best, effs):
    np.save('data/improved_angles', improved_angles)
    np.save('data/eff_best', eff_best)
    np.save('data/effs', np.array(effs))

def loop_through_cells(improved_angles, eff_best, effs, iter_num=None):
    for i in range(Nx):
        for j in range(Ny):
            if iter_num is not None:
                print(f'\niteration number = {iter_num+1}\n')
            improved_angles, eff_best = line_search(i, j, improved_angles, eff_best, effs)
    return improved_angles, eff_best

num_iterations = 3
improved_angles = starting_angles.copy()
eff_best = eff_design
effs = [eff_design]

# single pass through each angle
Nx, Ny = starting_angles.shape
print(f'will need to adjust {Nx * Ny} angles')
print(f'should take at least {(2*Nx*Ny*t_single*num_iterations/60/60):.2f} hours')

for iter_num in range(num_iterations):
    save_progress(improved_angles, eff_best, effs)
    improved_angles, eff_best = loop_through_cells(improved_angles, eff_best, effs, iter_num=iter_num)


f, (ax1, ax2, ax3) = plt.subplots(1, 3, tight_layout=True, figsize=(10, 5))
ax1.plot([100*e for e in effs])
ax1.set_xlabel('number of simulations')
ax1.set_ylabel('efficiency (%)')

ax2.plot(t_single / 60 / 60 * np.arange(1, len(effs)+1), [100*e for e in effs])
ax2.set_xlabel('total time (hours)')
ax2.set_ylabel('efficiency (%)')

t_solve = 4.06333  # from log

ax3.plot(t_solve  /60 * np.arange(1, len(effs)+1), [100*e for e in effs])
ax3.set_xlabel('solver time (min)')
ax3.set_ylabel('efficiency (%)')

plt.savefig('progress.png', dpi=400)

plot_angles(starting_angles)
plot_angles(improved_angles)

geometry = create_geometry(starting_angles)
sim_start = td.Simulation(size=sim_size,
                    mesh_step=[dl, dl, dl],
                    structures=geometry,
                    sources=source_circ,
                    monitors=monitors,
                    run_time=run_time,
                    pml_layers=pml_layers)
run(sim_start, task_name='spline_perturb')

geometry = create_geometry(improved_angles)
sim_end = td.Simulation(size=sim_size,
                    mesh_step=[dl, dl, dl],
                    structures=geometry,
                    sources=source_circ,
                    monitors=monitors,
                    run_time=run_time,
                    pml_layers=pml_layers)
run(sim_end, task_name='spline_perturb')

f, (ax1, ax2) = plt.subplots(1, 2, tight_layout=True)
sim_start.viz_field_2D(monitor_3fwhm, val='int', ax=ax1)
sim_end.viz_field_2D(monitor_3fwhm, val='int', ax=ax2)
plt.savefig('start_end.png', dpi=400)

