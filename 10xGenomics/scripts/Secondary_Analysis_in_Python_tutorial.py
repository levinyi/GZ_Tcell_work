'''This notebook demonstrates:

loading a feature-barcode matrix in HDF5 format
loading analysis files in CSV format
plotting UMI and feature (gene) count distributions
plotting clustering results and gene-specific expression in tSNE space
Data:

Go to the 10k Cells from a combined cortex, hippocampus and subventricular zone of an E18 mouse (v3 chemistry) page.
Download the Feature / cell matrix HDF5 (filtered) (.h5 format) and Clustering analysis (.tar format) files.
Move the downloaded files to the current foler.
Note:

This requires the following Python packages: Matplotlib, NumPy, SciPy, Pandas, h5py.
The easiest way to obtain these packages (as well as Jupyter Notebook) is to install Anaconda.'''

# import modules, define some functions for loading, saving and processing a gene-barcode matrix

import collections
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import h5py


