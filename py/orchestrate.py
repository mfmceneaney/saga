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

def get_list(divisions):
    """
    Parameters
    ----------
    divisions : dictionary, required
        Map of option names to option values

    Description
    -----------
    Get a list of all possible option value combinations from a map of option names to possible values.
    """

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = []
    for i, key in enumerate(divisions):
        if i==0:
                for el in divisions[key]:
                    data_list.append({key: el})
        else:
            data_list_new = []
            for d in data_list:
                for el in divisions[key]:
                    data_list_new.append(d.copy())
                    data_list_new[-1][key] = el
            data_list = data_list_new

    return data_list

def create_jobs(divisions,base_dir,submit_path,yaml_path):
    """
    Parameters
    ----------
    divisions : dictionary, required
        Map of option names to option values
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

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

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

def submit_jobs(divisions,base_dir,submit_path,out_path):
    """
    Parameters
    ----------
    divisions : dictionary, required
        Map of option names to option values
    base_dir : string, required
        Path to directory in which to create job directories
    submit_path : string, required
        Path to base version of SLURM job submission script
    yaml_path : string, required
        Path to base version of yaml file containing arguments for the executable run in the SLURM job submission script

    Description
    -----------
    Submit jobs within a job directory structure created by `orchesrate.create_jobs()`.
    """

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = get_list(divisions)

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
        command = command+'sbatch '+os.path.abspath(submit_path_i)+' >> '+out_path, #NOTE: COMMENTED OUT FOR DEBUGGING
        subprocess.run(command, shell=True, check=True, text=True)
        counter += 1
