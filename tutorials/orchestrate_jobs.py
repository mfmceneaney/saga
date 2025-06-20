from saga.orchestrate import create_jobs, submit_jobs

# Create job submission structure
methods = {"method":["HB","LF"]}
fitvars = {"asymfitvars":[["costheta1"],["costheta2"]]}
seeds   = {"inject_seed":[2**i for i in range(4)]}

# Results file paths and config
base_dir    = "results/"
submit_path = base_dir+"submit.sh"
yaml_path   = base_dir+"args.yaml"
out_path    = base_dir+"jobs.txt"
configs = dict(
    methods,
    **fitvars,
    **seeds
)
create_jobs(configs,base_dir,submit_path,yaml_path)
submit_jobs(configs,base_dir,submit_path,out_path,dry_run=True)

# Now write some dummy data for the purposes of the tutorial
import os
import numpy as np
from saga.aggregate import (
    get_config_list,
    get_config_str,
    get_out_file_name,
    get_out_dirs_list,
    get_binscheme_cuts_and_ids,
)
from saga.data import save_txt, load_yaml
def write_dummy_data_file_from_config(outdir):
    """
    WRITE DUMMY DATA FOR A CONFIGURATION ASSUMING THE :obj:`args.yaml` FILE EXISTS
    AND CONTAINS A BIN SCHEME NAMED :obj:`binscheme`.
    """

    binscheme_name = "binscheme" #NOTE: MAKE SURE THIS MATCHES THE DUMMY BIN SCHEME NAME IN ARGS.YAML

    yaml_path = os.path.join(outdir,"args.yaml")
    yaml_args = load_yaml(yaml_path)
    binscheme = yaml_args["binschemes"][binscheme_name]

    # Get the list of bin ids
    binscheme_cuts, binscheme_cut_titles, binscheme_ids, nested_grid_shape = get_binscheme_cuts_and_ids(
        binscheme,
        start_idx=0,
        id_key='bin_id'
    )

    # Get the name of the CSV file for the binning scheme you are interested in
    out_file_name = get_out_file_name(
            base_dir=outdir, #NOTE APPEND SYSTEM SEPARATOR TO OUTDIR HERE TODO: MAYBE JUST CHANGE THIS TO AN OUTPUT DIRECTORY...
            base_name='',
            binscheme_name=binscheme_name,
            ext='.csv'
        )

    # Convert to getKinBinnedAsym() output CSV format
    # COLS: bin_id,count,{binvarmean,binvarerr},{depolvarmean,depolvarerr},{asymfitvar,asymfitvarerr}

    # Set column header
    delimiter=","
    err_ext = "err"
    binvar_cols = [delimiter.join([binvar,binvar+err_ext]) for binvar in binscheme.keys()]
    depolvar_cols = [delimiter.join([depolvar,depolvar+err_ext]) for depolvar in ["depol"]]
    asymfitvar_cols = [delimiter.join([asymfitvar,asymfitvar+err_ext]) for asymfitvar in ["a0"]]
    cols = ["bin_id","count",*binvar_cols,*depolvar_cols,*asymfitvar_cols]
    header = delimiter.join(cols)

    # Set column formats
    fmt = ["%.3g" for i in range(2*(len(cols)-2))]
    fmt = ["%d","%d", *fmt]

    # Set bin means
    binvar_means = {binvar: np.array([np.average([binscheme[binvar][i],binscheme[binvar][i+1]]).item() for i in range(len(binscheme[binvar])-1)]) for binvar in binscheme}
    sum_binvar_nbins = np.sum([len(binvar_means[binvar]) for binvar in binscheme])

    # Set linearly dependent data with random noise
    data = []
    id_key = "bin_id"
    for binscheme_ids_idx in range(len(binscheme_ids)):

        # Set bin id and count
        bin_id = binscheme_ids[id_key].iloc[binscheme_ids_idx]
        count = 100

        # Set row data
        binvar_data = np.array([[binvar_means[binvar][binscheme_ids[binvar].iloc[binscheme_ids_idx]],0.05] for binvar in binscheme.keys()]).flatten().tolist()
        depolvar_data = [0.5,0.0]
        asymfitvar_val = -1.0 + 2.0 * np.sum([binscheme_ids[binvar].iloc[binscheme_ids_idx] for binvar in binscheme.keys()]) / sum_binvar_nbins  #NOTE: Introduce a linear dependence on each bin variable
        asymfitvar_val += np.random.rand() * 0.05 #NOTE: Introduce some random noise
        asymfitvar_data = [asymfitvar_val,0.05]
        row = [bin_id,count,*binvar_data,*depolvar_data,*asymfitvar_data]

        # Add row data
        data.append(row) 

    # Save data to file
    save_txt(
        out_file_name,
        data,
        delimiter=",",
        header=header,
        fmt=fmt,
        comments='',
    )

# Loop directories and write dummy data
out_dirs_list = get_out_dirs_list(configs,base_dir)
for out_dirs in out_dirs_list:
    for outdir in out_dirs:
        write_dummy_data_file_from_config(outdir)
