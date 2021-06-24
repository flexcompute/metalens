# metalens
Simulation and optimization of dielectric metalens with Tidy3D.

This repo contains an example of using Tidy3D to simulate and optimize a large area metalens.

![alt text](https://github.com/flexcompute/metalens/blob/main/Metalens.png?raw=true)

## Contents

### Simulation Files

- The jupyter notebook `Metalens_Simulate_Single.ipynb` runs the single simulation of a large area, 100 wavelength diameter metalens (with total dimensions of (100 x 100 x 46 wavelengths) in under 5 minutes.

- The python script `Metalens_Optimize.py` runs a brute force optimization of the metalens unit cells to improve focusing efficiency by 10%.


### Results and Post Processing

- The colormap for plotting has data stored in `BGYR_cmp.txt` and is imported for plotting from `cmap.py`.

- The data from the optimization runs are stored in `data` directory, as well as a script for potting this data.

## Finally

If you want more information on signing up to use tidy3d for free, please see our [documentation](https://simulation.cloud/docs/html/index.html) or reach out to us directly!
