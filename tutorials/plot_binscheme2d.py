import numpy as np
import os
import uproot as ur
import matplotlib.pyplot as plt
from matplotlib import colors

import saga.aggregate as sagas
from saga.data import load_yaml, load_th1
from saga.plot import (
    set_default_plt_settings,
    plot_th2,
    get_lims_coords,
    plot_lines,
    get_bin_centers,
    plot_bin_ids
)

# Setup, modify these as needed for your specific binning scheme
yaml_path = os.path.abspath('results_kinematics/args.yaml')
hist_path = os.path.abspath('out_binscheme_binvars_2D.root')
hist_name = 'h2__x_Q2'
binscheme_name = 'binscheme' #NOTE: Use `binscheme_grid` in the example yaml to plot some grid bin scheme limits.
binvars = ['x','Q2']
binvar_labels = {'x':'$x$','Q2':'$Q^{2}$ (GeV)$^{2}$'}
binvar_lims = {'x':[0.0,1.0],'Q2':[1.0,11.0]}
outpath = f'binscheme2d_{binvars[0]}_{binvars[1]}.pdf'
var_keys = []#binvars #NOTE: This should only be set in the case of a 2D grid scheme.
start_idx = 0
id_key = 'bin_id'

# Read bin scheme from YAML
yaml_args = load_yaml(yaml_path)
binscheme = yaml_args['binschemes'][binscheme_name]

# Load TH2 histogram with uproot
h2 = load_th1(hist_path,name=hist_name)

# Set plt settings
set_default_plt_settings()

# Open the figure
f, ax = plt.subplots(figsize=(16,10))

# Plot the 2D distribution
plot_th2(h2, ax, norm=colors.LogNorm())
ax.set_xlabel(binvar_labels[binvars[0]])
ax.set_ylabel(binvar_labels[binvars[1]])

# Get the bin limit line coordinates and plot
lims_coords = get_lims_coords(binscheme, binvar_lims['x'], binvar_lims['Q2'], var_keys=var_keys)
plot_lines(ax, lims_coords, linecolor='red', linewidth=1)

# Get bin scheme cuts and ids
cuts, _, _, _ = sagas.get_binscheme_cuts_and_ids(
                                                    binscheme,
                                                    start_idx=start_idx,
                                                    id_key=id_key,
                                                    binvar_titles=None,
                                                )

# Get the bin centers
bin_centers, bin_widths = get_bin_centers(cuts,swap_axes=False)

# Plot the bin ids
sagas.plot_bin_ids(
        ax,
        bin_centers,
        bin_widths=bin_widths,
        size=25,
        color='red',
        alpha=1.0,
    )

# Save the figure
f.savefig(outpath)
