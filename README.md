# metalens
Simulation and optimization of dielectric metalens with Tidy3D.

This repo contains the code and data for the paper `Simulation and Optimization of Large Area Metalenses`.

## Contents

### Simulation Files

- The jupyter notebook `Metalens_Simulate_Single.ipynb` runs the single simulation of the large area, 100 wavelength diameter metalens.

- The python script `Metalens_Optimize.py` runs the brute force optimization.

- Note that timings may vary between your results and the paper due to changes in the code over time.

### Results and Post Processing

- The colormap for plotting has data stored in `BGYR_cmp.txt` and is imported for plotting from `cmap.py`.

- The data from the optimization runs are stored in `data` directory, as well as a script for potting this data.

## Finally

- If you find this example useful, please cite our paper.

- If you want more information on signing up to use tidy3d for free, please see our [documentation](https://simulation.cloud/docs/html/index.html) or reach out to us directly!
