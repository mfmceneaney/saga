import numpy as np
import os
import pandas as pd

from saga.plot import get_bin_kinematics_title, plot_results_array

# Setup, modify these as needed for your specific binning scheme
csv_path = os.path.abspath('results_kinematics/out_binscheme_kinematics.csv')
hist_path = os.path.abspath('results_kinematics_hists/out_binscheme_kinematics.root')
grid_shape = (3,2)
kinvars = ['mass_pipim', 'phperp_pipim', 'z_pipim']
xlabels = {'mass_pipim':'$M_{\\pi^{+}\\pi^{-}}$ (GeV)','phperp_pipim':'$P_{\\perp, \\pi^{+}\\pi^{-}}$ (GeV)','z_pipim':'$z_{\\pi^{+}\\pi^{-}}$'}
xlims = {'mass_pipim':[0.0,3.0],'phperp_pipim':[0.0,1.25],'z_pipim':[0.0,1.0]}
hist_colors = {'mass_pipim':['tab:blue'],'phperp_pipim':['tab:red'],'z_pipim':['tab:orange']}

# Load kinematics CSV
df = pd.read_csv(csv_path)
bin_ids = df['bin'].unique().tolist()

# Set graph and plot_results arrays
graph_array = [[{} for j in range(len(bin_ids))] for i in range(len(kinvars))]
plot_results_kwargs_array = [[
        {
            'hist_keys':[f'h1_bin{bin_id}_'+kinvar],
            'title':get_bin_kinematics_title(bin_id,df),
            'xlims':xlims[kinvar],
            'xlabel':xlabels[kinvar],
            'hist_colors':hist_colors[kinvar],
        }
        for bin_id in bin_ids
    ]
    for kinvar in kinvars
]

# Set base kwargs
plot_results_kwargs_base = {
    'ylims':[0.0,450],
    'show_injected_asymmetries':False,
    'hist_clone_axis':False,
    'hist_paths':[hist_path],
    'hist_labels':['RGH MC'],
    'hist_linewidth':4,
    'ylabel': 'Counts',
    'watermark':'',
    'hist_density':False
}

# Set additional kwargs
figsize = (16*len(bin_ids),10*len(kinvars))
outpath = 'rgh_kinematics.pdf'
use_default_plt_settings = True
use_grid_xlabels = False

# Plot an array of graphs
plot_results_array(
        graph_array,
        plot_results_kwargs_array,
        plot_results_kwargs_base = plot_results_kwargs_base,
        figsize = figsize,
        outpath = outpath,
        use_default_plt_settings = use_default_plt_settings,
        use_grid_xlabels = use_grid_xlabels, #NOTE: Since you plot different x-axis variables in each row make sure that the labels are not dropped for rows below the top row.
    )
