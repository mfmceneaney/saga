import numpy as np
import os
import sys
import uproot as ur
import matplotlib.pyplot as plt
from matplotlib import colors

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../py')))

import saga.aggregate as sagas

# Setup, modify these as needed for your specific binning scheme
yaml_path = os.path.abspath('getBinKinematicsTH1Ds.yaml')
hist_path = os.path.abspath('out_binscheme_binvars_2D.root')
hist_name = 'h2__x_Q2'
binscheme_name = 'binscheme' #NOTE: Use `binscheme_grid` in the example yaml to plot some grid bin scheme limits.
binvars = ['x','Q2']
binvar_labels = {'x':'$x$','Q2':'$Q^{2}$ (GeV)$^{2}$'}
binvar_lims = {'x':[0.0,1.0],'Q2':[1.0,11.0]}
outpath = f'binscheme2d_{binvars[0]}_{binvars[1]}.pdf'
var_keys = []#binvars #NOTE: This should only be set in the case of a 2D grid scheme.

# Read bin scheme from YAML
yaml_args = sagas.load_yaml(yaml_path)
binscheme = yaml_args['binschemes'][binscheme_name]

# Load TH2 histogram with uproot
h2 = sagas.load_th1(hist_path,name=hist_name)

# Set plt settings
sagas.set_default_plt_settings()

# Open the figure
f, ax = plt.subplots(figsize=(16,10))

# Plot the 2D distribution
sagas.plot_TH2(h2, ax, norm=colors.LogNorm())
ax.set_xlabel(binvar_labels[binvars[0]])
ax.set_ylabel(binvar_labels[binvars[1]])

# Get the bin limit line coordinates and plot
lims_coords = sagas.get_lims_coords(binscheme, binvar_lims['x'], binvar_lims['Q2'], var_keys=var_keys)
sagas.plot_lines(ax, lims_coords, linecolor='red', linewidth=1)

# Save the figure
f.savefig(outpath)
