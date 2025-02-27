#----------------------------------------------------------------------#
# Python module for submitting slurm jobs for all possible
# option combinations supplied in yaml files.
# Authors: M. McEneaney (2024, Duke University)
#----------------------------------------------------------------------#
import subprocess
import os
import shutil
import yaml
import sys
import numpy as np
from .aggregate import get_config_list, get_config_str, get_binscheme_cuts_and_ids, save_txt, load_yaml, get_out_file_name

def write_file_from_config(outdir):
    """
    WRITE DUMMY DATA FOR A CONFIGURATION ASSUMING THE `args.yaml` FILE EXISTS
    AND CONTAINS A BIN SCHEME NAMED `binscheme.
    """

    binscheme_name = "binscheme" #NOTE: MAKE SURE THIS MATCHES THE DUMMY BIN SCHEME NAME IN ARGS.YAML

    yaml_path = os.path.join(outdir,"args.yaml")
    yaml_args = load_yaml(yaml_path)
    binscheme = yaml_args["binschemes"][binscheme_name]

    # Get the list of bin ids
    binscheme_cuts , binscheme_ids = get_binscheme_cuts_and_ids(
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

def create_jobs(configs,base_dir,submit_path,yaml_path):
    """
    Parameters
    ----------
    configs : dictionary, required
        Map of configuration option names to option values
    base_dir : string, required
        Path to directory in which to create job directories
    submit_path : string, required
        Path to base version of SLURM job submission script
    yaml_path : string, required
        Path to base version of yaml file containing arguments for the executable run in the SLURM job submission script

    Description
    -----------
    Create job directories for all possible option value combinations from a map of option names to possible values.
    Each directory will contain an appropriately modified version of the supplied SLURM and yaml files.
    """

    # Create map of elements of elements of configs and combine completely into each other for one list
    data_list = get_config_list(configs)

    # Loop resulting list
    for data_list_i in data_list:

        # Make job directory name and directory
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,"_".join([str(ele) for ele in data_list_i[key]]) if type(data_list_i[key])==list else str(data_list_i[key]) ]) for key in sorted(data_list_i)]))
        data_list_i["outdir"] = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
        try:
            os.mkdir(job_dir)
        except FileExistsError:
            print("WARNING: ",job_dir,"already exists.")

        # Copy files into job directory
        submit_path_i = os.path.join(job_dir,os.path.basename(submit_path))
        yaml_path_i   = os.path.join(job_dir,os.path.basename(yaml_path))
        cp_result = shutil.copyfile(submit_path,submit_path_i)
        shutil.copyfile(yaml_path,yaml_path_i)

        # Replace key values in yaml file from list and insert if non-existent
        with open(yaml_path_i, 'r') as yaml_i:
            doc = yaml.safe_load(yaml_i)
        for key in data_list_i:
            doc[key] = data_list_i[key]
        with open(yaml_path_i, 'w') as yaml_i:
            yaml.safe_dump(doc, yaml_i, default_flow_style=False)

        # Replace path to yaml with path to job yaml in submit script
        with open(submit_path_i, 'r') as submit_i:
            doc = submit_i.read()
        doc = doc.replace(os.path.abspath(base_dir),os.path.abspath(job_dir)) #NOTE: YAML SHOULD SPECIFY THE OUTPUT DIRECTORY FOR THE JOB. COULD HAVE A DEFAULT JOB DIRECTORY TO REPLACE THOUGH...
        with open(submit_path_i, 'w') as submit_i:
            submit_i.write(doc)

def submit_jobs(configs,base_dir,submit_path,out_path,dry_run=False):
    """
    Parameters
    ----------
    configs : dictionary, required
        Map of configuration option names to option values
    base_dir : string, required
        Path to directory in which to create job directories
    submit_path : string, required
        Path to base version of SLURM job submission script
    yaml_path : string, required
        Path to base version of yaml file containing arguments for the executable run in the SLURM job submission script
    dry_run : bool, optional
        Option to just print commands to text file and not submit jobs via sbatch
        Default : False

    Description
    -----------
    Submit jobs within a job directory structure created by `orchesrate.create_jobs()`.
    """

    # Create map of elements of elements of configs and combine completely into each other for one list
    data_list = get_config_list(configs)

    # Loop resulting list
    counter = 1
    for data_list_i in data_list:

        # Get job directory name
        job_dir = os.path.join(base_dir,"__".join(["_".join([key,"_".join([str(ele) for ele in data_list_i[key]]) if type(data_list_i[key])==list else str(data_list_i[key]) ]) for key in sorted(data_list_i)]))
        data_list_i["outdir"] = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.

        # Get submit script name
        submit_path_i = os.path.join(job_dir,os.path.basename(submit_path))
         
        # Submit job to SLURM
        command = 'echo \''+str(counter)+' sbatch '+os.path.abspath(submit_path_i)+'\' >> '+out_path+'; '
        if not dry_run: command = command+'sbatch '+os.path.abspath(submit_path_i)+' >> '+out_path
        subprocess.run(command, shell=True, check=True, text=True)
        if dry_run: write_file_from_config(job_dir)#NOTE: Write dummy data to configuration with a static binning scheme
        counter += 1
