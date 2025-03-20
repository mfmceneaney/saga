import numpy as np
import os
import sys
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../py')))

import saga.aggregate as sagas

# Setup, modify these as needed for your specific binning scheme
csv_path = os.path.abspath('results_kinematics_hists_2D/out_binscheme_kinematics.csv')
hist_path = os.path.abspath('results_kinematics_hists_2D/out_binscheme_kinematics_2d.root')
grid_shape = (3,2)
kinvars = [['z_pipim', 'phperp_pipim']]
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
            'hist_keys':[f'h2_bin{bin_id}_'+kinvar_x+'_'+kinvar_y],
            'title':sagas.get_bin_kinematics_title(bin_id,df),
            'xlims':xlims[kinvar_x],
            'xlabel':xlabels[kinvar_x],
            'ylims':xlims[kinvar_y],
            'ylabel':xlabels[kinvar_y],
        }
        for bin_id in bin_ids
    ]
    for kinvar_x, kinvar_y in kinvars
]

# Set base kwargs
plot_results_kwargs_base = {
    'show_injected_asymmetries':False,
    'hist_clone_axis':False,
    'hist_paths':[hist_path],
    'hist_labels':['RGH MC'],
    'watermark':'',
    'hist_density':False,
    'axlinewidth':0,
    'hist_dim':2,
    'legend_loc':None #NOTE: Do not plot a legend since you are using 2d hists.
}

# Set additional kwargs
figsize = (16*len(bin_ids),10*len(kinvars))
outpath = 'rgh_kinematics_2d.pdf'
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
