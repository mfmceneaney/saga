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

# Set up chaining for batched data (specifically `old_dat_path`)
nbatch = 1
nbatches = {"nbatches":[nbatch]}
ibatches = {"ibatch":[i for i in range(nbatch)]}
chain_keys = ["nbatches", "ibatch"]
chain_configs = dict(
    nbatches,
    **ibatches,
) if nbatch > 1 else {}
aggregate_config = {"inject_seed":1} if nbatch > 1 else {} #NOTE: You must set this to correctly determine the path when chaining and aggregating.

# Setup input paths
base_dir     = os.path.abspath("results/")
submit_path  = os.path.join(base_dir,"submit.sh")
yaml_path    = os.path.join(base_dir,"args.yaml")
out_path     = os.path.join(base_dir,"jobs.txt")

# Set aggregate keys
aggregate_keys = ["inject_seed"]

# Load the binscheme you want to use
binschemes_name = "binschemes"
binscheme_name = 'binscheme'
yaml_args = sagas.load_yaml(yaml_path)
binscheme = yaml_args[binschemes_name][binscheme_name]

# Arguments for sagas.get_config_list()
result_name = "a0" #NOTE: This also gets recycled as the asymmetry name

# Arguments for sagas.get_out_dirs_list()
sep='_'
ext='.pdf'

# Arguments for sagas.get_out_file_name()
out_file_name_ext = '.csv'

# Arguments for sagas.rescale_csv_data()
rescale_csv_data_kwargs_base = {
        'old_dat_path':'results/',
        'new_sim_path':'results/',
        'old_sim_path':'results/',
        'count_key':'count',
        'yerr_key':result_name,
        'yerr_key':result_name+'err',
        'xs_ratio':1.0,
        'lumi_ratio':1.0,
        'tpol_factor':0.5,
        'tdil_factor':1.0,
        'yvalue':-100.0,
}

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
    out_dirs = out_dirs_list[config_idx]

    # Get the name of the CSV file for the binning scheme you are interested in
    out_file_names = [sagas.get_out_file_name(
            base_dir=outdir,
            binscheme_name=binscheme_name,
            ext=out_file_name_ext
        ) for outdir in out_dirs]

    # Rescale results files
    for out_file_name in out_file_names:
        sagas.rescale_csv_data(
            out_file_name,
            outpath = '', #NOTE: Output path will be first argument with `_rescaled` inserted before the file extension if this is empty
            config=config_list[config_idx],
            aggregate_config=aggregate_config,
            chain_configs=chain_configs,
            **rescale_csv_data_kwargs_base,
        )
