import numpy as np
import matplotlib
paper_cmap_colors = np.loadtxt('BGYR_cmp.txt', delimiter=',')
paper_cmap = matplotlib.colors.ListedColormap(paper_cmap_colors)