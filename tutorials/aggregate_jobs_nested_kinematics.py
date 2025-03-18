import numpy as np
import os
import sys
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../py')))

import saga.aggregate as sagas

# Setup, modify these as needed for your specific binning scheme
grid_shape = (3,2)
kinvars = ['mass_pipim', 'phperp_pipim', 'z_pipim']
xlabels = {'mass_pipim':'$M_{\pi^{+}\pi^{-}}$ (GeV)','phperp_pipim':'$P_{\perp, \pi^{+}\pi^{-}}$ (GeV)','z_pipim':'$z_{\pi^{+}\pi^{-}}$'}
xlims = {'mass_pipim':[0.0,3.0],'phperp_pipim':[0.0,1.25],'z_pipim':[0.0,1.0]}
hist_colors = {'mass_pipim':['tab:blue'],'phperp_pipim':['tab:red'],'z_pipim':['tab:orange']}

# Load kinematics CSV
csv_path = "/path/to/out_binscheme_kinematics.csv" #NOTE: Modify this to the appropriate path
df = pd.read_csv(csv_path)
bin_ids = df['bin'].unique().tolist()

# Set graph and plot_results arrays
graph_array = [[{},{},{},{},{}] for i in range(len(kinvars))]
plot_results_kwargs_array = [[
        {
            'hist_keys':[f'h1_bin{bin_id}_'+kinvar],
            'title':sagas.get_bin_kinematics_title(bin_id,df),
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
    'ylims':[0.0,800],
    'show_injected_asymmetries':False,
    'hist_clone_axis':False,
    'hist_paths':['/path/to/out_binscheme_kinematics.root'], #NOTE: Modify this to the appropriate path
    'hist_labels':['RGH MC'],
    'ylabel': 'Counts',
    'watermark':'',
    'hist_density':False
}

# Set additional kwargs
figsize = (16*grid_shape[0]*grid_shape[1],10*len(kinvars))
outpath = 'rgh_kinematics.pdf'
use_default_plt_settings = True
use_grid_xlabels = False

# Plot an array of graphs
sagas.plot_results_array(
        graph_array,
        plot_results_kwargs_array,
        plot_results_kwargs_base = plot_results_kwargs_base,
        figsize = figsize,
        outpath = outpath,
        use_default_plt_settings = use_default_plt_settings,
        use_grid_xlabels = use_grid_xlabels, #NOTE: Since you plot different x-axis variables in each row make sure that the labels are not dropped for rows below the top row.
    )
