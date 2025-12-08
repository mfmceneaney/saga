"""
This module implements methods for submitting
slurm jobs for all possible option combinations
supplied in yaml argument files.
"""

import subprocess
import os
import shutil
import yaml
from .aggregate import get_config_list, get_config_str


def create_jobs(
    configs, base_dir="./", submit_path="submit.sh", yaml_path="args.yaml", aliases=None, replacements=None,
):
    """
    Parameters
    ----------
    configs : dict, required
        Map of configuration option names to option values
    base_dir : str, required
        Path to directory in which to create job directories
    submit_path : str, required
        Path to base version of SLURM job submission script
    yaml_path : str, required
        Path to base version of yaml file containing arguments for
        the executable run in the SLURM job submission script
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases
    replacements : dict, optional
        Map of yaml keys to default values, overriding those in :obj:`yaml_path` and :obj:`configs`

    Description
    -----------
    Create job directories for all possible option value
    combinations from a map of option names to possible values.
    Each directory will contain an appropriately modified version
    of the supplied SLURM and yaml files.
    """

    # Get config list
    data_list = get_config_list(configs)

    # Loop resulting list
    for data_list_i in data_list:

        # Make job directory name and directory
        job_dir = os.path.join(base_dir, get_config_str(data_list_i, aliases=aliases))
        data_list_i["outdir"] = os.path.abspath(job_dir)
        try:
            os.mkdir(job_dir)
        except FileExistsError:
            print("WARNING: ", job_dir, "already exists.")

        # Copy files into job directory
        submit_path_i = os.path.join(job_dir, os.path.basename(submit_path))
        yaml_path_i = os.path.join(job_dir, os.path.basename(yaml_path))
        shutil.copyfile(submit_path, submit_path_i)
        shutil.copyfile(yaml_path, yaml_path_i)

        # Replace key values in yaml file from list and insert if non-existent
        with open(yaml_path_i, "r", encoding="utf-8") as yaml_i:
            doc = yaml.safe_load(yaml_i)
        for key in data_list_i:
            doc[key] = data_list_i[key]
        if isinstance(replacements, dict):
            for key in replacements:
                doc[key] = replacements[key]
        with open(yaml_path_i, "w", encoding="utf-8") as yaml_i:
            yaml.safe_dump(doc, yaml_i, default_flow_style=False)

        # Replace path to yaml with path to job yaml in submit script
        with open(submit_path_i, "r", encoding="utf-8") as submit_i:
            doc = submit_i.read()
        doc = doc.replace(os.path.abspath(base_dir), os.path.abspath(job_dir))
        with open(submit_path_i, "w", encoding="utf-8") as submit_i:
            submit_i.write(doc)


def submit_jobs(
    configs,
    base_dir="./",
    submit_path="submit.sh",
    out_path="jobs.txt",
    aliases=None,
    dry_run=False,
):
    """
    Parameters
    ----------
    configs : dict, required
        Map of configuration option names to option values
    base_dir : str, required
        Path to directory in which to create job directories
    submit_path : str, required
        Path to base version of SLURM job submission script
    yaml_path : str, required
        Path to base version of yaml file containing arguments for
        the executable run in the SLURM job submission script
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases
    dry_run : bool, optional
        Option to just print commands to text file and not submit jobs via sbatch

    Description
    -----------
    Submit jobs within a job directory structure created by :meth:`orchestrate.create_jobs`.
    """

    # Get config list
    data_list = get_config_list(configs)

    # Loop resulting list
    counter = 1
    for data_list_i in data_list:

        # Get job directory name
        job_dir = os.path.join(base_dir, get_config_str(data_list_i, aliases=aliases))
        data_list_i["outdir"] = os.path.abspath(job_dir)

        # Get submit script name
        submit_path_i = os.path.join(job_dir, os.path.basename(submit_path))

        # Submit job to SLURM
        command = "echo '" + str(counter)
        command += " sbatch " + os.path.abspath(submit_path_i)
        command += "' >> " + out_path + "; "
        if not dry_run:
            command = (
                command + "sbatch " + os.path.abspath(submit_path_i) + " >> " + out_path
            )
        subprocess.run(command, shell=True, check=True, text=True)
        counter += 1
