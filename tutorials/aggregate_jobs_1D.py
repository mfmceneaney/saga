import numpy as np
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'../py')))

import saga.aggregate as sagas
import matplotlib.pyplot as plt

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
base_dir     = os.path.abspath("results_1D/")
submit_path  = os.path.join(base_dir,"submit.sh")
yaml_path    = os.path.join(base_dir,"args.yaml")
out_path     = os.path.join(base_dir,"jobs.txt")
#NOTE: Set the bin migration path below since this is binscheme dependent

# Set aggregate keys
aggregate_keys = ["inject_seed"]

# Load the binschemes
binschemes_name = "binschemes"
yaml_args = sagas.load_yaml(yaml_path)
binschemes = yaml_args[binschemes_name]

# Arguments for sagas.get_config_list()
result_name = "a0" #NOTE: This also gets recycled as the asymmetry name

# Arguments for sagas.get_config_out_path()
sep='_'
ext='.pdf'

# Arguments for sagas.get_out_file_name()
out_file_name_ext = '.csv'
bin_mig_base_name="bin_mig_mat_"

# Arguments for sagas.apply_bin_mig()
use_bin_mig = True
id_gen_key='binid_gen'
id_rec_key='binid_rec'
mig_key='mig'
results_keys = [result_name] #NOTE: You can apply bin migration to multiple dataframe entries in one go.

# Arguments for sagas.get_graph_data()
id_key = 'bin_id'

# Arguments for sagas.get_graph_data()
count_key  = 'count'
asym_key   = result_name #NOTE: This is set from above
err_ext    = 'err'

# Arguments for sagas.plot_results()
plot_results_kwargs_base = {
    'ylims':[-1.0,2.0],
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

# Additional useful parameters for plotting
xlabel_map = {binscheme_name:'x' for binscheme_name in binschemes}#NOTE: Just map bin scheme names to x variable titles for ease of use
hist_paths_map = {binscheme_name:['results/h1.root'] for binscheme_name in binschemes}
figsize = (16,10)
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
            'xs_ratio':1.0,
            'lumi_ratio':1.0,
        },
    )

#---------- Set configurations ----------#
# Get list of configurations
config_list = sagas.get_config_list(configs,aggregate_keys=aggregate_keys)

# Get aggregated list of directories
out_dirs_list = sagas.get_out_dirs_list(
                                configs,
                                base_dir,
                                aggregate_keys=aggregate_keys
                            )

#---------- Loop bin schemes ----------#
for binscheme_idx, binscheme_name in enumerate(binschemes.keys()):

    # Get the bin scheme
    binscheme = binschemes[binscheme_name]
    proj_var  = list(binscheme.keys())[0] #NOTE: Assume projection variable is the only variable in the bin scheme
    nbins = len(binscheme[proj_var])-1

    # Arguments for sagas.get_graph_data()
    xvar_keys = [proj_var]

    # Set some bin scheme dependent plotting parameters
    binlims = binscheme[proj_var]
    plot_results_kwargs_base['binlims'] = binlims
    plot_results_kwargs_base['xlims'] = [binlims[0],binlims[-1]]
    plot_results_kwargs_base['xlabel'] = xlabel_map[binscheme_name]
    plot_results_kwargs_base['hist_paths'] = hist_paths_map[binscheme_name]

    # Get the bin migration path
    bin_mig_path = sagas.get_out_file_name(
        base_dir=base_dir,
        base_name=bin_mig_base_name,
        binscheme_name=binscheme_name,
        ext=out_file_name_ext
    )

    # Load bin migration matrix and invert
    bin_mig_df, bin_mig_mat, inv_bin_mig_mat = None, None, None
    if use_bin_mig:
        bin_mig_df = sagas.load_csv(bin_mig_path)
        bin_mig_mat = sagas.get_bin_mig_mat(
            bin_mig_df,
            id_gen_key=id_gen_key,
            id_rec_key=id_rec_key,
            mig_key=mig_key,
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

    #---------- Loop configurations ----------#
    # Loop each aggregate list
    for config_idx in range(len(config_list)):

        # Set the config you are interested in
        config = config_list[config_idx]
        out_dirs = out_dirs_list[config_idx]

        # Set the output path basename for this config
        config_out_path = sagas.get_config_out_path(
                base_dir,
                aggregate_keys,
                binscheme_name+sep+result_name,
                config,
                sep=sep,
                ext=ext,
            )
        config_out_path = os.path.join(base_dir,config_out_path)

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

        # Get an aggregate graph
        proj_ids = [i for i in range(nbins)]#NOTE: Assume bin scheme indices are simple
        sgasym_idx = plot_results_kwargs_base['sgasym_idx'] #NOTE: Assume this is in the kwargs base dictionary
        sgasym = config['sgasyms'][sgasym_idx] if 'sgasyms' in config else 0.0
        aggregate_graph = sagas.get_aggregate_graph(
            [
                sagas.get_graph_data(
                            df,
                            proj_ids,
                            id_key=id_key,
                            count_key=count_key,
                            xvar_keys=xvar_keys,
                            asym_key=asym_key,
                            err_ext=err_ext
                ) for df in dfs
            ],
            xvar_keys=xvar_keys,
            sgasym=sgasym
        )

        # Use default plotting settings
        if use_default_plt_settings: sagas.set_default_plt_settings()

        # Create figure and axes
        f, ax = plt.subplots(figsize=figsize)

        # Set additional arguments for saga.aggregate.plot_results()
        plot_results_kwargs_base['sgasyms'] = config['sgasyms']
        plot_results_kwargs_base['outpath'] = config_out_path

        # Plot the graph
        sagas.plot_results(ax,**aggregate_graph,**plot_results_kwargs_base)

        # Save the graph
        f.savefig(config_out_path)
