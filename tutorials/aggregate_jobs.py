import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../py')))

import saga.aggregate as sagas

# Setup configuration dictionary
methods = {"method":["HB","LF"]}
fitvars = {"asymfitvars":["costheta1","costheta2"]}
seeds   = {"inject_seed":[2**i for i in range(4)]}
configs = dict(
    methods,
    **fitvars,
    **seeds,
)

# Setup input paths
base_dir     = os.path.abspath("results/")
submit_path  = os.path.join(base_dir,"submit.sh")
yaml_path    = os.path.join(base_dir,"args.yaml")
out_path     = os.path.join(base_dir,"jobs.txt")
bin_mig_path = os.path.join(base_dir,"bin_mig_mat.csv")

# Set aggregate keys
aggregate_keys = ["inject_seed"]

# Load the binscheme you want to use
binschemes_name = "binschemes"
binscheme_name = 'binscheme'
yaml_args = sagas.load_yaml(yaml_path)
binscheme = yaml_args[binschemes_name][binscheme_name]

# Load bin migration matrix and invert
use_bin_mig = True
id_gen_key='binid_gen'
id_rec_key='binid_rec'
mig_key='mig'
bin_mig_df, bin_mig_mat, inv_bin_mig_mat = None, None, None
if use_bin_mig:
    bin_mig_df = sagas.load_csv(bin_mig_path)
    bin_mig_mat = sagas.get_bin_mig_mat(
        bin_mig_df,
        id_gen_key='binid_gen',
        id_rec_key='binid_rec',
        mig_key='mig',
    )
    sagas.save_bin_mig_mat_to_csv(
        bin_mig_mat,
        base_dir='./',
        basename=binscheme_name,
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    )
    inv_bin_mig_mat = np.linalg.inv(bin_mig_mat)

# Arguments for sagas.get_config_list()
result_name = "a0" #NOTE: This also gets recycled as the asymmetry name

# Arguments for sagas.get_out_dirs_list()
sep='_'
ext='.pdf'

# Arguments for sagas.get_out_file_name()
out_file_name_ext = '.csv'

# Arguments for sagas.apply_bin_mig()
results_keys = [result_name] #NOTE: You can apply bin migration to multiple dataframe entries in one go.

# Arguments for sagas.get_binscheme_cuts_and_ids()
id_key = 'bin_id'
start_idx = 0
binvar_titles = ['x','y','z']

# Arguments for sagas.get_projection_ids()
proj_vars = ['x']
arr_vars = ['y','z']
arr_var_bins={} #NOTE: Only set this if you want to restrict yourself to specific bins in `arr_vars`, but you will need to correctly reset the `grid_shape` below.

# Arguments for sagas.get_graph_data()
count_key  = 'count'
asym_key   = result_name #NOTE: This is set from above
err_ext    = 'err'
xvar_keys  = proj_vars

# Arguments for sagas.get_graph_array()
sgasym = 0.1

# Arguments for sagas.plot_results_array()
proj_vars_shape = [len(binscheme[proj_var])-1 for proj_var in proj_vars]
arr_vars_shape = [len(binscheme[arr_var])-1 for arr_var in arr_vars]
grid_shape = [*arr_vars_shape,*proj_vars_shape]
plot_results_kwargs_array = [[
    {
        'ylims':[-1.0,2.0],
        'xlims':[0.0,4.0],#NOTE: You can set kwargs depending on the location in the grid here
    } for j in range(grid_shape[1])] for i in range(grid_shape[0])
]
binlims = binscheme[xvar_keys[0]]
plot_results_kwargs_base = {
    'binlims':binlims,
    'sgasyms':[0.1,0.05,-0.3],
    'sgasym_idx':0,
    'sgasym_labels':['$\mathcal{A}_{1}$','$\mathcal{A}_{2}$','$\mathcal{A}_{3}$'],
    'sg_colors':['blue','red','green'],
    'bgasyms':[0.1],
    'bgasym_labels':['$\mathcal{A}$'],
    'bg_colors':['cyan'],
    'show_injected_asymmetries':True,
    'hist_paths':['results/h1.root'],
    'hist_colors':['tab:orange'],
    'hist_keys':['h1'],
    'hist_labels':['MC distribution'],
}
figsize = (16*grid_shape[0],10*grid_shape[1])
#NOTE: Set outpath within the loop for unique naming
use_default_plt_settings = True

# If you want to rescale your results using results from other base directories set the following arguments
rescale = False
if rescale:
    plot_results_kwargs_base = dict(
        plot_results_kwargs_base,
        **{
            'old_dat_path':'results/',
            'new_sim_path':'results/',
            'old_sim_path':'results/',
            'count_key':'count',
            'yerr_key':'',
            'lumi_ratio':1.0,
        },
    )

#---------- Set configurations and loop ----------#

# Get list of configurations
config_list = sagas.get_config_list(configs,aggregate_keys=aggregate_keys)

# Get aggregated list of directories
out_dirs_list = sagas.get_out_dirs_list(
                                configs,
                                base_dir,
                                aggregate_keys=aggregate_keys
                            )

# Loop each aggregate list
for config_idx in range(len(config_list)):

    # Set the config you are interested in
    config = config_list[config_idx]
    out_dirs = out_dirs_list[config_idx]

    # Set the output path basename for this config
    config_out_path = sagas.get_config_out_path(
            base_dir,
            aggregate_keys,
            result_name,
            config,
            sep=sep,
            ext=ext,
        )

    # Get the name of the CSV file for the binning scheme you are interested in
    out_file_names = [sagas.get_out_file_name(
            base_dir=outdir,
            binscheme_name=binscheme_name,
            ext=out_file_name_ext
        ) for outdir in out_dirs]

    # Load pandas dataframes from the files
    dfs = [sagas.load_csv(out_file_name) for out_file_name in out_file_names]

    # Apply bin migration correction
    if use_bin_mig:
        for df in dfs:
            sagas.apply_bin_mig(df,inv_bin_mig_mat,results_keys=results_keys) #NOTE: THIS MODIFIES THE DATAFRAMES IN PLACE

    # Get bin scheme cuts and ids
    binscheme_cuts, binscheme_cut_titles, binscheme_ids = sagas.get_binscheme_cuts_and_ids(
                                                        binscheme,
                                                        start_idx=start_idx,
                                                        id_key=id_key,
                                                        binvar_titles=binvar_titles,
                                                    )

    # Get projection bin ids
    all_proj_ids, arr_vars, all_proj_arr_var_ids = sagas.get_projection_ids(
            binscheme_ids,
            proj_vars,
            arr_vars = arr_vars,
            id_key=id_key,
            arr_var_bins=arr_var_bins
        )

    # Open a single graph
    graph_data = sagas.get_graph_data(
                            dfs[0],
                            all_proj_ids[0][0],
                            id_key=id_key,
                            count_key=count_key,
                            xvar_keys=xvar_keys,
                            asym_key=asym_key,
                            err_ext=err_ext
                )

    # Open an array of graphs with the shape of the projection ids array
    graph_array = sagas.get_graph_array(
            dfs,
            all_proj_ids,
            id_key=id_key,
            count_key=count_key,
            xvar_keys=xvar_keys,
            asym_key=asym_key,
            err_ext=err_ext,
            sgasym=sgasym,
        )

    # Get array of bin cut titles
    cut_array = sagas.get_cut_array(
        binscheme_cut_titles,
        all_proj_ids,
        arr_vars,
    )

    # Modify `plot_results_kwargs_array` setting bin cuts as titles
    sagas.add_cut_array(
        plot_results_kwargs_array,
        cut_array,
        arr_vars,
    )

    # Plot an array of graphs
    sagas.plot_results_array(
            graph_array,
            plot_results_kwargs_array,
            plot_results_kwargs_base = plot_results_kwargs_base,
            figsize = figsize,
            outpath = config_out_path,
            use_default_plt_settings = use_default_plt_settings,
        )
