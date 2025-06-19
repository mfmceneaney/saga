#----------------------------------------------------------------------#
# Python module for aggregating output from slurm jobs for all possible
# option combinations supplied in yaml files.
# Authors: M. McEneaney (2024, Duke University)
#----------------------------------------------------------------------#
import uproot as ur
import numpy as np
import pandas as pd
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sbn

import subprocess
import os
import shutil
import yaml
import sys

def get_config_str(
        config,
        sep='_',
        aliases={}
    ):
    """
    Parameters
    ----------
    config : dict, required
        Configuration mapping option names to option values
    sep : str, optional
        String separator
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases

    Returns
    -------
    str
        String representation of configuration

    Description
    -----------
    Create a string representation of a configuration.  Note that aliases of nonhashable types, e.g., :obj:`list` or :obj:`dict`
    will be accessed by the string representation of the aliased object :obj:`str(<object>)`.
    """

    return (sep+sep).join([
                aliases[key][config[key]] if (key in aliases and (type(config[key]) in (str,float,int)) and config[key] in aliases[key])
                else aliases[key][str(config[key])] if (key in aliases and (type(config[key]) not in (str,float,int)) and str(config[key]) in aliases[key])
                else sep.join([
                    key,sep.join([str(ele) for ele in config[key]]) if type(config[key])==list else str(config[key])
                ]) for key in sorted(config)
                ])

def get_config_out_path(
        base_dir,
        aggregate_keys,
        result_name,
        config,
        sep='_',
        aliases={},
        ext='.pdf',
    ):
    """
    Parameters
    ----------
    base_dir : str, required
        Directory name to prepend to output file names
    aggregate_keys : list, required
        List of keys over which configurations are grouped
    result_name : str, required
        Unique string identifier for result
    config : dict, required
        Job configuration map
    sep : str, optional
        String separator
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases
    ext : str, optional
        File extension

    Returns
    -------
    str
        Unique output path for a file produced from the given configuration

    Description
    -----------
    Get a unique output file name for a given configuration.
    """

    job_config_name  = 'aggregate'+sep+sep+sep+(sep+sep).join([str(key) for key in sorted(aggregate_keys)])+(sep+sep+sep)
    job_config_name += get_config_str(config,sep=sep,aliases=aliases)
    job_config_name += (sep+sep+sep)+result_name+ext
    outpath = os.path.abspath(os.path.join(base_dir,job_config_name))

    return outpath

def get_out_file_name(
        base_dir = None,
        base_name='',
        binscheme_name='x',
        ext='.csv'
    ):
    """
    Parameters
    ----------
    base_dir : str, optional
        Base directory for file path
    base_name : str, optional
        Base name used to construct file name
    binscheme_name : str, optional
        Name of binning scheme
    ext : str, optional
        File extension

    Returns
    -------
    str
        Output file name

    Description
    -----------
    Get output file name for a generic kinematic binning scheme passed to an executable constructed as :obj:`baseoutpath+binscheme_name+ext`.
    """
    out_file_name = ''.join([base_name,binscheme_name,ext])
    return out_file_name if base_dir is None else os.path.join(base_dir,out_file_name)

def get_config_list(
        configs,
        aggregate_keys=[]
    ):
    """
    Parameters
    ----------
    configs : dict, required
        Map of configuration option names to option values
    aggregate_keys : list, optional
        List of keys over which to group configurations

    Description
    -----------
    Get a list of all possible option value combinations from a map of option names to possible values.
    """

    # Create map of elements of elements of configs and combine completely into each other for one list
    data_list = []
    first = True
    for key in configs:
        if key in aggregate_keys: continue
        if first:
            for el in configs[key]:
                data_list.append({key: el})
            first = False
        else:
            data_list_new = []
            for d in data_list:
                for el in configs[key]:
                    data_list_new.append(d.copy())
                    data_list_new[-1][key] = el
            data_list = data_list_new

    return data_list

def get_out_dirs_list(
        configs,
        base_dir,
        aggregate_keys=[],
        aliases={}
    ):
    """
    Parameters
    ----------
    configs : dict, required
        Map of configuration option names to option values
    base_dir : str, required
        Path to directory in which to create job directories
    aggregate_keys : list, optional
        List of keys over which to group configurations
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases

    Returns
    -------
    list
        List of job output directory names grouped across the values for :obj:`aggregate_keys`

    Description
    -----------
    Get a list of output directory names grouped across the keys listed in :obj:`aggregate_keys` for jobs generated
    from all possible option value combinations from a map of configuration option names to possible values.
    """
    
    # Get a list of all possible option value combinations from configs
    config_list = get_config_list(configs,aggregate_keys=aggregate_keys)
    
    # Create list for output directories lists
    out_dirs_list = []

    # Loop resulting list
    for config_list_i in config_list:

        # Create output directories list
        output_dirs = []

        # Check if aggregate keys is empty
        if len(aggregate_keys)==0:

            # Get job directory
            config_str = get_config_str(config_list_i,aliases=aliases)
            job_dir = os.path.abspath(os.path.join(base_dir,config_str))
            output_dirs.append(job_dir)
            
        # Loop aggregate keys and build file list for current binning
        for key in aggregate_keys:
            for value in configs[key]:
            
                # Set the aggregate key value
                config_list_i_val = dict(config_list_i) #NOTE: Clone dictionary so that you only edit the local one.
                config_list_i_val[key] = value

                # Get job directory
                config_str = get_config_str(config_list_i_val,aliases=aliases)
                job_dir = os.path.abspath(os.path.join(base_dir,config_str))
                output_dirs.append(job_dir)

        # Now add the list of output directories to the overall list
        out_dirs_list.append(output_dirs)

    return out_dirs_list

def load_th1(
        path,
        name = "h1"
    ):
    """
    Parameters
    ----------
    path : str, required
        Path to ROOT file containing histogram
    name : str, optional
        Name of :obj:`TH1` object within the ROOT file

    Returns
    -------
    np.array
        Histogram data as a numpy array or empty list if file is not found

    Description
    -----------
    Read :obj:`TH1` histogram data from a ROOT file.  Note that this will work for any histogram: (:obj:`TH1`, :obj:`TH2`, :obj:`TH3`).
    """

    # Get TH1 from ROOT file
    if path=="": return []
    try:
        f = ur.open(path)
        g = f[name].to_numpy()
        return g

    except FileNotFoundError as e:
        print("FileNotFoundError: ",path)
        print("\t Returning empty list")
        return []

def load_yaml(
        path
    ):
    """
    Parameters
    ----------
    path : str, required
        Path to yaml file

    Returns
    -------
    dict
        dictionary of yaml contents

    Description
    -----------
    Load a yaml from a file.
    """

    yaml_args = {}
    with open(path) as f:
        yaml_args = yaml.safe_load(f)
    return yaml_args

def load_csv(
        path,
        old_path=None,
        new_path=None,
        config={},
        aggregate_config={},
        chain_configs={},
        aliases={},
    ):
    """
    Parameters
    ----------
    path : str, required
        Path to CSV file containing table
    old_path : str, optional
        Directory name to replace
    new_path : str, optional
        Directory name to insert if not None
    config : dict, optional
        Map of configuration option names to option values
    aggregate_config : dict, optional
        Map of aggregate configuration option names to option values, used for determining correct directory names
    chain_configs : dict, optional
        Map of configuration option names to lists of values across which to chain
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe of csv file data

    Description
    -----------
    Read table stored in CSV format into a pandas DataFrame optionally replacing
    part of the path, for example, a directory name, with another value.
    """
    inpath = path.replace(old_path,new_path) if old_path is not None and new_path is not None else path

    # Chain CSVs across the given keys
    if len(config)>0 and len(chain_configs)>0:

        # Get the full batch config
        aggregate_configs = {key:[aggregate_config[key]] for key in aggregate_config}
        configs = dict(
            {key:[config[key]] for key in config},
            **chain_configs,
            **aggregate_configs
        )

        # Get a list of all possible option value combinations from configs
        config_list = get_config_list(configs,aggregate_keys=[])

        # Set csv list
        csv_list = []

        # Loop resulting list
        for config_list_i in config_list:

            # Get job directory
            config_str = get_config_str(config_list_i,aliases=aliases)

            # Get base job directory
            base_config_str = get_config_str(config,aliases=aliases)

            # Modify path for chain element
            inpath_i = inpath.replace(base_config_str, config_str)

            # Open csv
            csv_i = pd.read_csv(inpath_i)
            csv_list.append(csv_i)

        # Merge csvs and return
        return pd.concat(csv_list, ignore_index=True)

    # Return csv
    else: return pd.read_csv(inpath)

def set_nested_bin_cuts(
        cuts,
        cut_titles,
        ids,
        node,
        old_cuts       = [],
        old_cut_titles = [],
        old_ids        = [],
        var            = "",
        lims_key       = "lims",
        nested_key     = "nested",
        binvar_titles  = {}
    ):

    """
    Parameters
    ----------
    cuts : list, required
        List of cuts to recursively update
    cut_titles : list, required
        List of cut titles to recursively update
    ids : list, required
        List of bin ids to recursively update
    old_cuts : list, optional
        List of cuts from previous recursion
    old_cut_titles : list, optional
        List of cut titles from previous recursion
    old_ids : list, optional
        List of bin ids from previous recursion
    var : str, optional
        Bin scheme variable at current recursion depth
    lims_key : str, optional
        Key for bin limits
    nested_key : str, optional
        Key for nested bin scheme
    binvar_titles : dict, optional
        Map of bin scheme variable names to LaTeX titles

    Description
    -----------
    Recursively set lists of bin cuts, titles, and ids for each bin scheme variable
    in a nested binning scheme.
    """

    # Check the YAML node
    if (type(node)==dict):

        # Check for bin limits
        lims = []
        if (lims_key in node.keys() and type(node[lims_key])==list):
            lims = node[lims_key]

        # Set nbins lower limit to 0 since you allow passing limits with length 0
        nbins = len(lims)-1
        if (nbins<0): nbins=0

        # Get bin variable title
        var_title = var if var not in binvar_titles else binvar_titles[var]

        # Loop bins and get bin cuts
        varlims = [[lims[idx],lims[idx+1]] for idx in range(nbins)]
        varids  = [i for i in range(nbins)]
        varcuts = [f"({var}>={varlims[idx][0]} && {var}<{varlims[idx][1]})" for idx in range(nbins)]
        varcut_titles = [f"${varlims[idx][0]} \\leq {var_title} < {varlims[idx][1]}$" for idx in range(nbins)]

        # Expand bin cuts and cut titles and projection ids maps
        old_cuts = [f"{varcut} && {cut}" for varcut in varcuts for cut in old_cuts] if len(old_cuts)>0 else varcuts
        old_cut_titles = [f"{var} : {varcut_title} && {cut_title}" for varcut_title in varcut_titles for cut_title in old_cut_titles] if len(old_cut_titles)>0 else [f"{var} : {varcut_title}" for varcut_title in varcut_titles]
        old_ids = [dict({var:varid},**id_map) for varid in varids for id_map in old_ids] if len(old_ids)>0 else [{var:el} for el in varids]

        # Check for nested binning
        if (nested_key in node.keys() and type(node[nested_key])==list):

            # Get nested YAML node
            node_nested = node[nested_key]

            # Loop nested bins
            for ibin in range(len(node_nested)): #NOTE: These are not bin limits just maps to bin limits for each bin, so loop normally.

                # Get nested YAML node
                node_nested_bin = node_nested[ibin]

                # Loop nested bin variables (only expect one!)
                for it_key in node_nested_bin:

                    # Get nested YAML node
                    it_nested = node_nested_bin[it_key]

                    # Create a new vector for uniqueness along different recursion branches
                    new_old_cuts = []
                    new_old_cut_titles = []
                    new_old_ids = []
                    if (ibin<len(old_cuts)):
                        new_old_cuts = [old_cuts[ibin]]
                        new_old_cut_titles = [old_cut_titles[ibin]]
                        new_old_ids = [old_ids[ibin]]

                    # Recursion call
                    set_nested_bin_cuts(
                        cuts,
                        cut_titles,
                        ids,
                        it_nested,
                        new_old_cuts,
                        new_old_cut_titles,
                        new_old_ids,
                        it_key,
                        lims_key,
                        nested_key
                    )

                    # Break on first nested variable found
                    break

        # Case you don't have further nested binning
        else:
            for ibin in range(len(old_cuts)):
                cuts.append(old_cuts[ibin])
                cut_titles.append(old_cut_titles[ibin])
                ids.append(old_ids[ibin])
            return
    else:
        return

def get_scheme_vars(
        node,
        nested_key = "nested"
    ):
    """
    Parameters
    ----------
    node : dict, required
        Map of bin scheme variables to bin limits arrays with either a nested or grid structure
    nested_key : str, optional
        Key for nested bin scheme structure

    Description
    -----------
    Find the bin scheme variables in a bin scheme with either a nested or grid structure.
    """

    # Initialize array
    binscheme_vars = []

    # Check for nested bin scheme
    if (type(node)==dict and nested_key in node.keys() and type(node[nested_key])==list and type(node[nested_key][0])==dict):

        # Get nested node
        node_nested = dict(node)

        while (type(node_nested)==dict and nested_key in node_nested.keys() and type(node_nested[nested_key])==list and type(node_nested[nested_key][0])==dict):

            binvar = list(node_nested[nested_key][0].keys())[0]
            binscheme_vars.append(binvar)
            node_nested = dict(node_nested[nested_key][0][binvar])
        
        return binscheme_vars

    # Default to case you have a grid scheme
    else:
        return [var for var in node.keys()]

def get_nested_scheme_shape(
        node,
        lims_key   = "lims",
        nested_key = "nested"
    ):
    """
    Parameters
    ----------
    node : dict, required
        Map of bin scheme variables to bin limits arrays with a nested structure
    lims_key : str, optional
        Key for bin limits arrays
    nested_key : str, optional
        Key for nested bin scheme structure

    Raises
    ------
    ValueErrror:
        Raise an error if :obj:`node` does not define a 2D nested bin structure.

    Returns
    -------
    list
        List of lengths of nested bin variable bins for each bin in the outer variable
        of a 2D nested bin scheme

    Description
    -----------
    Find the bin scheme shapes for a 2D bin scheme with nested structure.
    """

    # Initialize array
    binscheme_shape = []

    # Check for nested bin scheme
    if (type(node)==dict and nested_key in node.keys() and type(node[nested_key])==list and type(node[nested_key][0])==dict):

        # Get nested node
        node_nested = dict(node)

        # Loop nested node
        depth = 0
        while (type(node_nested)==dict and nested_key in node_nested.keys() and type(node_nested[nested_key])==list and type(node_nested[nested_key][0])==dict):
            
            # Get bin variable
            binvar = list(node_nested[nested_key][0].keys())[0]
            
            # Increment and cheeck depth
            depth += 1
            if depth>=2:

                # Add shapes to list at a recursion depth of 2 and exit loop
                for idx in range(len(node_nested[nested_key])):
                    shape = 0
                    if lims_key in node_nested[nested_key][idx][binvar].keys():
                        shape = len(node_nested[nested_key][idx][binvar][lims_key])-1
                        shape = max(0,shape)
                    binscheme_shape.append(shape)
                break
            
            # Update nested node
            node_nested = dict(node_nested[nested_key][0][binvar])
        
        return binscheme_shape

    # Case you do not have a nested bin scheme
    else:
        raise ValueError("`node` must define have a 2D nested bin scheme structure")

def get_binscheme_cuts_and_ids(
        binscheme,
        start_idx=0,
        id_key='bin_id',
        binvar_titles = [],
    ):
    """
    Parameters
    ----------
    binscheme : dict, required
        Dictionary mapping binning variables to bin limits arrays
    start_idx : int, optional
        Starting integer for enumerating bin unique integer ids
    id_key : str, optional
        Column name for bin unique integer ids
    binvar_titles : list, optional
        List of latex format titles for bin variables for forming cut titles

    Returns
    -------
    (dict, dict, pandas.DataFrame, list or None)
        Dictionary of bin ids to cuts, dictionary of bin ids to bin cut
        title dictionaries, a pandas dataframe with bin ids under
        :obj:`id_key` and the projection bin ids under the respective variable names,
        and the shape of the nested bin scheme at depth 2.  Note that the nested
        grid shape will be :obj:`None` if no 2D nested bin scheme is defined.

    Description
    -----------
    Create a maps of bin ids to bin cuts and titles and a data frame mapping bin ids to projection bin ids.
    Also, check and return the nested bin scheme shapes at depth 2 in the case of a 2D nested bin scheme.
    """

    # Initialize arrays
    cuts = []
    cut_titles = []
    ids  = []

    # Initialize nested shape
    nested_grid_shape = None

    # Set titles using raw bin variable names if no titles are provided
    binscheme_vars = get_scheme_vars(binscheme,nested_key='nested')
    if binvar_titles is None or len(binvar_titles)!=len(binscheme_vars): binvar_titles = [var for var in binscheme_vars]

    # Check for nested bin scheme
    if (type(binscheme)==dict and 'nested' in binscheme.keys() and type(binscheme['nested'])==list and type(binscheme['nested'][0])==dict):

        # Recursively set bin scheme cuts map
        set_nested_bin_cuts(
            cuts,
            cut_titles,
            ids,
            binscheme,
            lims_key      = 'lims',
            nested_key    = 'nested',
            binvar_titles = binvar_titles
        )

        nested_grid_shape = get_nested_scheme_shape(
            binscheme,
            lims_key   = 'lims',
            nested_key = 'nested'
        )

    # Grid bin scheme case
    else:

        # Loop bin variables
        for var_idx, var in enumerate(binscheme):

            # Set the variable title
            var_title  = binvar_titles[var_idx]

            # Get bin variable limits, projection ids, cuts and cut titles
            nbins = len(binscheme[var])-1
            lims = binscheme[var]
            if type(binscheme[var])==dict and 'nbins' in binscheme[var].keys() and 'lims' in binscheme[var].keys():
                nbins = binscheme[var]['nbins']
                lims = [(binscheme[var]['lims'][1]-binscheme[var]['lims'][0])/nbins * i + binscheme[var]['lims'][0] for i in range(nbins+1)]
            varlims = [[lims[idx],lims[idx+1]] for idx in range(nbins)]
            varids  = [i for i in range(nbins)]
            varcuts = [f"({var}>={varlims[idx][0]} && {var}<{varlims[idx][1]})" for idx in range(nbins)]
            varcut_titles = [f"${varlims[idx][0]} \\leq {var_title} < {varlims[idx][1]}$" for idx in range(nbins)]

            # Expand bin cuts and cut titles and projection ids maps
            cuts = [f"{varcut} && {cut}" for varcut in varcuts for cut in cuts] if len(cuts)>0 else varcuts
            cut_titles = [f"{var} : {varcut_title} && {cut_title}" for varcut_title in varcut_titles for cut_title in cut_titles] if len(cut_titles)>0 else [f"{var} : {varcut_title}" for varcut_title in varcut_titles]
            ids  = [dict({var:varid},**id_map) for varid in varids for id_map in ids] if len(ids)>0 else [{var:el} for el in varids]

    # Turn cuts and cuts_titles into maps
    cuts = {start_idx+idx: cuts[idx] for idx in range(len(cuts))}
    cut_titles = {start_idx+idx: {el.split(" : ")[0]:el.split(" : ")[1] for el in cut_titles[idx].split(" && ")} for idx in range(len(cut_titles))}

    # Set up data frame
    df = {id_key:[]}
    for var in binscheme_vars:
        df[var] = []

    # Add data frame entries
    for idx in range(len(ids)):
        binscheme_binid = idx + start_idx
        df[id_key].append(binscheme_binid)
        for var in binscheme_vars:
            df[var].append(ids[idx][var])

    #  Create pandas data frame
    df = pd.DataFrame(df)

    return cuts, cut_titles, df, nested_grid_shape

def reshape_nested_grid(
        proj_ids,
        nested_grid_shape,
        fill_value = None,
    ):
    """
    Parameters
    ----------
    grid : list, required
        List of projection ids to reshape
    nested_grid_shape : list, required
        List of (irregular) nested grid array dimensions
    fill_value : int, optional
        Fill value for padding

    Returns
    -------
    list
        A reshaped grid array padded to dimension :obj:`(len(nested_grid_shape),max(nested_grid_shape))`

    Raises
    ------
    ValueError
        If :obj:`nested_grid_shape` is empty or not a list of integers or if :obj:`sum(nested_grid_shape)`
        is not the same length as the number of projections in :obj:`proj_ids`

    Description
    -----------
    Reshape projection ids for an irregular nested bin scheme padding to the largest
    nested dimension with :obj:`fill_value`.
    """

    # Get the shape of the new array and instantiate it
    shape = [len(nested_grid_shape),max(nested_grid_shape)]
    reshaped_proj_ids = [[fill_value for i in range(shape[1])] for j in range(shape[0])]

    # Check the nested grid shape
    if nested_grid_shape is None:
        return proj_ids
    if type(nested_grid_shape)!=list:
        raise TypeError('`nested_grid_shape` must be a list of integers')
    if len(nested_grid_shape)==0:
        raise ValueError('`nested_grid_shape` must be a non-empty list')
    if not type(nested_grid_shape[0])==int:
        raise TypeError('`nested_grid_shape` must be a list of integers')
    if sum(nested_grid_shape)!=len(proj_ids):
        raise ValueError('`sum(nested_grid_shape)` must be the same length as the number of projections in `proj_ids`')

    # Loop projection ids and set your new grid array values
    counter = 0
    for i, dim in enumerate(nested_grid_shape):
        for j in range(dim):
            reshaped_proj_ids[i][j] = proj_ids[counter]
            counter += 1

    return reshaped_proj_ids

def get_projection_ids(
        df,
        proj_vars,
        arr_vars = [],
        id_key='bin_id',
        arr_var_bins={},
        nested_grid_shape=None,
    ):
    """
    Parameters
    ----------
    df : pandas.DataFrame, required
        Pandas dataframe with bin ids under :obj:`id_key` and the projection bin ids under the respective variable names
    proj_vars : list, required
        Projection variables
    arr_vars : list, optional
        Variables in which to construct a grid of results
    id_key : str, optional
        Column name for bin unique integer ids
    arr_var_bins : dict, optional
        Dictionary of array binning scheme variable names to the specific projection bin ids to look for in those variables
    nested_grid_shape : list, optional
        One dimensional list of nested grid dimensions, note this will be padded with :obj:`None` to the largest nested dimension

    Returns
    -------
    (list, list, list)
        An array of each projection's unique bin ids, a list of the array binning variables encountered,
        and an array of the array variables' bin ids for each projection.  Note the array shapes should be :obj:`(*N_BINS_PROJ_VARS,*N_BINS_ARR_VARS)`.
        If :obj:`N_BINS_PROJ_VARS=[5]` and :obj:`N_BINS_ARR_VARS=[3,8]` then the shape of :obj:`all_proj_ids` and :obj:`all_proj_arr_var_ids`
        will  be :obj:`(3,8,5)`.

    Raises
    ------
    TypeError
        Raise an error if array binning variables found in keys of :obj:`arr_var_bins` are also
        found in :obj:`proj_vars` since should not simultaneously project and select a single bin.

    Description
    -----------
    Create an array of unique bin ids projected over a subset of binning variables from a binning scheme
    and organized in array-like structure over another subset of binning variables.
    """

    # Check projection variables
    if len(proj_vars)>2:
        print('WARNING: `get_projection_ids` : Are you sure you want more than 2 projection variables?')

    # Get list of grouping variables
    for key in df.keys():
        if (key!=id_key and key not in proj_vars and key not in arr_var_bins.keys()) and len(arr_vars)==0: arr_vars.append(key)

    # Check array variable bins argument
    if np.any([el in proj_vars for el in arr_vars]):
        raise TypeError('`proj_vars` entries are not allowed in `arr_vars`.')

    # Check array variable bins argument
    if np.any([el in proj_vars for el in arr_var_bins]):
        raise TypeError('`proj_vars` entries are not allowed in `arr_var_bins`.')

    # Get cut for array variable bins
    arr_var_bins_cut = ' and '.join([f'{key}=={arr_var_bins[key]}' for key in arr_var_bins])

    # Get unique projection bin ids for each projection variable
    unique_bin_ids = [df[var].unique() for var in proj_vars]
    nbins = [len(el) for el in unique_bin_ids]
    unique_arr_bin_ids = [df[var].unique() for var in arr_vars]
    nbins_arr = [len(el) for el in unique_arr_bin_ids]
    grid_shape = [*nbins_arr,*nbins]

    # Get map of starting bins and bin and projection indices
    query = ' and '.join([f'{proj_var}=={0}' for proj_var in proj_vars])
    if arr_var_bins_cut!='': query = ' and '.join([query,arr_var_bins_cut])
    start_bins_slice = df.query(query)

    # Loop starting bins and create projection lists
    all_proj_ids = []
    all_proj_arr_var_ids = []
    for proj_idx in range(len(start_bins_slice)):

        # Find bins that match starting bin in array bin variable indices
        arr_var_ids = [start_bins_slice[arr_var].iloc[proj_idx] for arr_var in arr_vars]
        query = ' and '.join([f'{arr_vars[idx]}=={arr_var_ids[idx]}' for idx in range(len(arr_var_ids))])
        if arr_var_bins_cut!='': query = ' and '.join([query,arr_var_bins_cut])
        proj_ids = df.query(query)
        proj_ids = proj_ids.sort_values(proj_vars)[id_key].values #NOTE: SORT BY PROJECTION VARIABLE BINS
        proj_ids = proj_ids.reshape(nbins) #NOTE: RESIZE BY APPROPRIATE NUMBER OF BINS IF REQUESTED.
        all_proj_ids.append(proj_ids.tolist())
        all_proj_arr_var_ids.append(arr_var_ids)

    all_proj_ids = np.reshape(all_proj_ids,grid_shape) if nested_grid_shape is None else reshape_nested(all_proj_ids, nested_grid_shape)
    all_proj_arr_var_ids = np.reshape(all_proj_arr_var_ids,(*nbins_arr,len(arr_vars)))

    return all_proj_ids, arr_vars, all_proj_arr_var_ids

def get_graph_data(
                    df,
                    bin_ids,
                    id_key='bin_id',
                    count_key='count',
                    xvar_keys=['x'],
                    asym_key='a0',
                    err_ext='err'
    ):
    """
    Parameters
    ----------
    df : pandas.DataFrame, required
        Pandas dataframe containing bin ids and data
    bin_ids : list, required
        List of unique integer bin ids
    id_key : str, optional
        String identifier for bin id column
    count_key : str, optional
        String identifier for bin count column
    xvar_keys : list, optional
        List of binning variables for which to return mean values
    asym_key : str, optional
        Asymmetry variables for which to return mean value
    err_ext : str, optional
        Extension for forming error column names

    Returns
    -------
    np.array
        Numpy array containing graph data with dimensions :obj:`(1+2*(1+len(xvar_keys)),*shape(bin_ids))`

    Description
    -----------
    Read graph data :obj:`(count, y, yerr, x0, x0err, x1, x1_err,...)` for a projection plot from a pandas dataframe.
    """

    # Initialize arrays
    cts  = []
    y    = []
    yerr = []
    x    = []
    xerr = []

    # Check if bin_ids is multi-dimensional
    bin_ids_shape = np.shape(bin_ids)
    if len(bin_ids_shape)>1:
        bin_ids = np.array(bin_ids).flatten()

    # Loop bins
    for bin_raw_idx, bin_id in enumerate(bin_ids):

        # Get bin data
        bin_data = df.loc[df[id_key]==bin_id]

        # Get bin count
        cts.append(bin_data[count_key].item())

        # Get bin asymmetry value and error
        y.append(bin_data[asym_key].item())
        yerr.append(bin_data[asym_key+err_ext].item())

        # Loop bin variables
        for xvar_idx, xvar_key in enumerate(xvar_keys):

            # Get bin variable data
            bin_x = bin_data[xvar_key].item()
            bin_x_err = bin_data[xvar_key+err_ext].item()#NOTE< MAAYBE USSE .LOC HERE???!?!?!?

            # Add bin variable mean and error
            if bin_raw_idx==0:
                x.append([bin_x])
                xerr.append([bin_x_err])
            else:
                x[xvar_idx].append(bin_x)
                xerr[xvar_idx].append(bin_x_err)

    # Reshape data
    if len(bin_ids_shape)>1:

        # Reshape bin counts
        cts   = np.reshape(cts,bin_ids_shape)

        # Reshape bin asymmetry value and error
        y    = np.reshape(y,bin_ids_shape)
        yerr = np.reshape(yerr,bin_ids_shape)

        # Reshape bin variable statistics
        for xvar_idx, xvar_key in enumerate(xvar_keys):
            x[xvar_idx]    = np.reshape(x[xvar_idx],bin_ids_shape)
            xerr[xvar_idx] = np.reshape(xerr[xvar_idx],bin_ids_shape)

    return np.array([
        cts,
        y,
        yerr,
        *x,
        *xerr
    ])

def get_aggregate_graph(
        graph_list,
        xvar_keys=['x'],
        sgasym=0.0,
    ):
    """
    Parameters
    ----------
    graph_list : list, required
        List of graphs with dimension :obj:`(N_GRAPHS, 2*(1+N_XVAR_KEYS), N_BIN_IDS)`
    xvar_keys : list, optional
        List of binning variables for which to return mean values
    sgasym : float, optional
        Injected signal asymmetry for computing difference of measured and injected values

    Returns
    -------
    dict
        Dictionary of mean asymmetry means and errors and other statistics names as well as bin variable means and errors to an array of their values in each kinematic bin

    Description
    -----------
    Compute the mean bin variables and asymmetry means and errors and other statistical information across a list of graphs' data from :meth:`get_graph_data`.
    Note that in the case of that the graph dimension is greater than 1, the bin variable statistics will be returned as a list of arrays in the same order as :obj:`xvar_keys`.
    """

    # Setup return dictionary
    graph = {
            'ct_mean':[],
            'y_mean':[],
            'yerr_mean':[],
            'y_std':[],
            'y_min':[],
            'y_max':[],
            'ydiff_mean':[],
            'ydiff_std':[],
            'ydiff_min':[],
            'ydiff_max':[],
            'x_mean':[],
            'xerr_mean':[]
            }

    # Check if graph list is empty
    if len(graph_list)==0:
        print("WARNING: len(graph_list)==0.  Returning empty graph.")
        return graph

    # Format graph list
    graph_list  = np.array(graph_list)
    graph_list  = np.swapaxes(graph_list,0,1) #NOTE: INCOMING LIST SHOULD GO FROM DIMENSION (N_GRAPHS, 2*(1+N_XVAR_KEYS), N_BIN_IDS) -> (2*(1+N_XVAR_KEYS), N_GRAPHS, N_BIN_IDS)

    # Extract aggregate asymmetry statistics
    ct_idx   = 0
    y_idx    = 1
    yerr_idx = 2
    graph['ct_mean']      = np.mean(graph_list[ct_idx],axis=0)
    graph['y_mean']     = np.mean(graph_list[y_idx],axis=0)
    graph['yerr_mean']  = np.sqrt(np.mean(np.square(graph_list[yerr_idx]),axis=0))
    graph['y_std']      = np.std(graph_list[y_idx],axis=0) if len(graph_list[y_idx])>1 else None
    graph['y_min']      = np.min(graph_list[y_idx],axis=0) if len(graph_list[y_idx])>1 else None
    graph['y_max']      = np.max(graph_list[y_idx],axis=0) if len(graph_list[y_idx])>1 else None
    graph['ydiff_mean'] = np.mean(graph_list[y_idx]-sgasym,axis=0)
    graph['ydiff_std']  = np.std(graph_list[y_idx]-sgasym,axis=0)
    graph['ydiff_min']  = np.min(graph_list[y_idx]-sgasym,axis=0)
    graph['ydiff_max']  = np.max(graph_list[y_idx]-sgasym,axis=0)

    # Extract aggregate projection variable statistics
    x_idx_start    = 3
    xerr_idx_start = x_idx_start+len(xvar_keys)
    if len(xvar_keys)==1:

        # Set projection variable arrays in the case of a 1D binning
        x_idx = x_idx_start
        xerr_idx = xerr_idx_start
        graph['x_mean']     = np.mean(graph_list[x_idx],axis=0) #NOTE: Get mean across different graphs (axis=0) but not across bins (axis>0)
        graph['xerr_mean']  = np.sqrt(np.mean(np.square(graph_list[xerr_idx]),axis=0))
    else:

        # Loop projection variable keys and aggregate across each variable for a >1D binning
        for xvar_idx, xvar in enumerate(xvar_keys):
            x_idx = x_idx_start + xvar_idx
            xerr_idx = xerr_idx_start + xvar_idx
            graph['x_mean'].append(np.mean(graph_list[x_idx],axis=0))
            graph['xerr_mean'].append(np.sqrt(np.mean(np.square(graph_list[xerr_idx]),axis=0)))

    return graph

def get_graph_array(
        dfs,
        proj_ids,
        id_key='bin_id',
        count_key='count',
        xvar_keys=['x'],
        asym_key='a0',
        err_ext='err',
        sgasym=0.0,
    ):
    """
    Parameters
    ----------
    dfs : list, required
        Array of pandas dataframes containing bin ids and data
    proj_ids : list, required
        Array of unique integer bin ids, note if entries are set to None the graph array entry will also be set to None to allow for masked grids
    id_key : str, optional
        String identifier for bin id column
    count_key : str, optional
        String identifier for bin count column
    xvar_keys : list, optional
        List of binning variables for which to return mean values
    asym_key : str, optional
        Asymmetry variables for which to return mean value
    err_ext : str, optional
        Extension for forming error column names
    sgasym : float or list, optional
        Injected signal asymmetry for computing difference of measured and injected values

    Returns
    -------
    list
        An array of aggregated graph data dictionaries obtained with the same grid structure as :obj:`proj_ids`.

    Raises
    ------
    TypeError
        Raise an error the shape of the array of bin indices is not in :math:`(1,2)`.

    Description
    -----------
    Aggregate a set of graphs given the array of dataframes and a grid array of projection bin indices.
    Note that the returned array will have the same shape as the given :obj:`proj_ids`.
    """

    # Get grid shape from projection ids array
    shape = np.shape(proj_ids)

    # Create a graph array in the 2D grid case
    if len(shape)==3:
        return [[
            get_aggregate_graph(
                [
                    get_graph_data(
                                df,
                                proj_ids[i][j],
                                count_key=count_key,
                                xvar_keys=xvar_keys,
                                asym_key=asym_key,
                                err_ext=err_ext
                    ) for df in dfs
                ],
                xvar_keys=xvar_keys,
                sgasym=sgasyms[i][j] if type(sgasym) is not float else sgasym
            ) if proj_ids[i][j] is not None else None for j in range(shape[1])] for i in range(shape[0]) #NOTE Allow masked grid
        ]

    # Create a graph array in the 1D grid case
    elif len(shape)==2:
        return [
            get_aggregate_graph(
                [
                    get_graph_data(
                                df,
                                proj_ids[i],
                                id_key=id_key,
                                count_key=count_key,
                                xvar_keys=xvar_keys,
                                asym_key=asym_key,
                                err_ext=err_ext
                    ) for df in dfs
                ],
                xvar_keys=xvar_keys,
                sgasym=sgasym[i] if type(sgasym) is not float else sgasym
            ) for i in range(shape[0])
        ]

    # Raise an error if another shape length is encountered
    else:
        raise TypeError('`get_graph_array` : `proj_ids` must have len(shape) in (2,3) but shape = ',shape)

def rescale_graph_data(
        ct_mean,
        x_mean,
        y_mean,
        xerr_mean,
        yerr_mean,
        path,
        old_dat_path = 'old_dat_path.csv',
        new_sim_path = 'new_sim_path.csv',
        old_sim_path = 'old_sim_path.csv',
        count_key = 'count',
        yerr_key = '',
        xs_ratio = 1.0,
        lumi_ratio = 1.0,
        tpol_factor = 1.0,
        tdil_factor = 1.0,
        yvalue = -100.0,
        xvar_keys = ['x'],
        sgasym = 0.0,
        aliases={},
    ):
    """
    Parameters
    ----------
    ct_mean : list, required
        Count mean values for each bin
    x_mean : list, required
        x mean values for each bin
    y_mean : list, required
        y mean values for each bin
    xerr_mean : list, required
        x error alues for each bin
    yerr_mean : list, required
        y error values for each bin
    path : str, required
        Path where the given graph data will be stored
    old_dat_path : str, optional
        Part of :obj:`path` to replace
    new_sim_path : str, optional
        Path with which to replace :obj:`old_dat_path` to find new simulation graph CSVs
    old_sim_path : str, optional
        Path with which to replace :obj:`old_dat_path` to find new simulation graph CSVs
    count_key : str, optional
        CSV column key for graph bin counts
    yerr_key : str, optional
        CSV column key for graph bin y errors
    xs_ratio : str, optional
        Cross-section ratio (new/old) for scaling
    lumi_ratio : str, optional
        Luminosity ratio (new/old) for scaling
    tpol_factor : float, optional
        Target polarization factor for rescaling
    tdil_factor : float, optional
        Target dilution factor for rescaling
    yvalue : float, optional
        Constant asymmetry value to be used for computing rescaled errors
    xvar_keys : list, optional
        List of binning variables for which to return mean values
    sgasym : float or list, optional
        Injected signal asymmetry for computing difference of measured and injected values
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases

    Returns
    -------
    dict
        A dictionary of graph data produced by :meth:`get_aggregate_graph`

    Description
    -----------
    Rescale a graph from data (:obj:`old_dat`) loading new and old simulation graphs (:obj:`new_sim` and :obj:`old_sim`)
    from file to compute the bin dependent acceptance ratio.  Start from the ratio (:obj:`new_sim`/:obj:`old_sim`) of counts in each bin.
    The rescaling ratio is then computed by multiplying by :obj:`lumi_ratio / xs_ratio`.  This ratio is used to rescale the counts, but
    the errors are further multiplied by a factor :math:`\\frac{1}{P_{Target} \\cdot D_{Target}}` to account for the target polarization
    and dilution factors.  Note that the asymmetry values can be set to a constant with :obj:`yvalue=A` if you only care about the rescaled errors,
    and that if you do so, the rescaled asymmetry errors will be further scaled as
    :math:`\\sigma_{A} = \\sqrt{\\frac{1-(A \\cdot P_{Target})^{2}}{N_{Rescaled}}}`.
    """

    # Load other graphs from csv
    new_sim_graph = load_csv(path,old_path=old_dat_path,new_path=new_sim_path,aliases=aliases)
    old_sim_graph = load_csv(path,old_path=old_dat_path,new_path=old_sim_path,aliases=aliases)

    # Get counts OR y errors from csv
    new_sim_graph_count = new_sim_graph[count_key] if yerr_key is None or yerr_key == '' else 1.0/np.square(new_sim_graph[yerr_key])
    old_sim_graph_count = old_sim_graph[count_key] if yerr_key is None or yerr_key == '' else 1.0/np.square(old_sim_graph[yerr_key])

    # Compute scaled quantities
    acceptanceratio  = np.divide(new_sim_graph_count,old_sim_graph_count) / xs_ratio
    acceptanceratio[np.isinf(acceptanceratio)] = 0
    acceptanceratio[np.isnan(acceptanceratio)] = 0 #NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    scaling          = acceptanceratio * lumi_ratio
    scaled_ct_mean   = np.multiply(scaling,ct_mean)
    err_scaling      = np.sqrt(np.divide(ct_mean,scaled_ct_mean)) #NOTE: SCALE ERRORS ASSUMING POISSONIAN STATISTICS -> d ~ 1/sqrt(N) -> MULTIPLY BY sqrt(N_old_data/N_new_data)
    err_scaling[np.isinf(err_scaling)] = 0
    err_scaling[np.isnan(err_scaling)] = 0 #NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    scaled_yerr_mean = np.multiply(err_scaling,yerr_mean) * 1.0/(tpol_factor * tdil_factor)

    # Set other graph quantities to new_sim values
    x_mean = new_sim_graph['x']
    xerr_mean = new_sim_graph['xerr']
    y_mean = new_sim_graph['y']

    # Set y values to constant and update scaled y errors if requested
    scaled_y_mean = y_mean if yvalue<-1 else [yvalue for i in range(len(y_mean))]
    if yvalue>=-1: scaled_yerr_mean *= np.sqrt(1-np.square(yvalue*tpol_factor))

    # Create a length 1 list of graph data with scaled graph results
    graph_list = np.array([[scaled_ct_mean,scaled_y_mean,scaled_yerr_mean,x_mean,xerr_mean]])

    graph = get_aggregate_graph(
            graph_list,
            xvar_keys=xvar_keys,
            sgasym=sgasym
        )

    # Set systematic errors to scaling fractions
    graph['scaling'] = scaling
    graph['acceptanceratio'] = acceptanceratio

    return graph

def rescale_csv_data(
        path,
        outpath = '',
        old_dat_path = 'old_dat_path.csv',
        new_sim_path = 'new_sim_path.csv',
        old_sim_path = 'old_sim_path.csv',
        count_key = 'count',
        y_key = 'a0',
        yerr_key = 'a0err',
        xs_ratio = 1.0,
        lumi_ratio = 1.0,
        tpol_factor = 1.0,
        tdil_factor = 1.0,
        yvalue = -100.0,
        float_format = "%.3g",
        config = {},
        aggregate_config = {},
        chain_configs = {},
        aliases={},
    ):
    """
    Parameters
    ----------
    path : str, required
        Path where the given graph data will be stored
    outfile_name : str, optional
        Output path, if empty the :obj:`path` argument will be used with :obj:`_rescaled` inserted before the :obj:`.csv` extension
    old_dat_path : str, optional
        Part of :obj:`path` to replace
    new_sim_path : str, optional
        Path with which to replace :obj:`old_dat_path` to find new simulation graph CSVs
    old_sim_path : str, optional
        Path with which to replace :obj:`old_dat_path` to find new simulation graph CSVs
    count_key : str, optional
        CSV column key for graph bin counts
    y_key : str, optional
        CSV column key for the asymmetry value corresponding to the rescaled errors
    yerr_key : str, optional
        CSV column key for graph bin y errors
    xs_ratio : str, optional
        Cross-section ratio (new/old) for scaling
    lumi_ratio : str, optional
        Luminosity ratio (new/old) for scaling
    tpol_factor : float, optional
        Target polarization factor for rescaling
    tdil_factor : float, optional
        Target dilution factor for rescaling
    yvalue : float, optional
        Constant asymmetry value to be used for computing rescaled errors
    float_format : str or Callable, optional
        Format string for floating point numbers passed to :meth:`pd.DataFrame.to_csv()`
    config : dict, optional
        Map of configuration option names to option values for chaining across :obj:`old_dat_path` CSVs
    aggregate_config : dict, optional
        Map of aggregate configuration option names to option values for chaining across :obj:`old_dat_path` CSVs
    chain_configs : dict, optional
        Map of configuration option names to lists of values across which to chain for :obj:`old_dat_path` CSVs
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases

    Description
    -----------
    Rescale a set of results from data (:obj:`old_dat`) loading new and old simulation results (:obj:`new_sim` and :obj:`old_sim`)
    from file to compute the bin dependent acceptance ratio.  Start from the ratio (:obj:`new_sim`/:obj:`old_sim`) of counts in each bin.
    The rescaling ratio is then computed by multiplying by :obj:`lumi_ratio / xs_ratio`.  This ratio is used to rescale the counts, but
    the errors are further multiplied by a factor :math:`\\frac{1}{P_{Target} \\cdot D_{Target}}` to account for the target polarization
    and dilution factors.  Note that the asymmetry values can be set to a constant with :obj:`yvalue=A` if you only care about the rescaled errors,
    and that if you do so, the rescaled asymmetry errors will be further scaled as
    :math:`\\sigma_{A} = \\sqrt{\\frac{1-(A \\cdot P_{Target})^{2}}{N_{Rescaled}}}`.
    """

    # Load results from csv
    old_dat_df = load_csv(path,config=config,aggregate_config=aggregate_config,chain_configs=chain_configs,aliases=aliases)
    new_sim_df = load_csv(path,old_path=old_dat_path,new_path=new_sim_path,aliases=aliases)#TODO: Could add other arguments for chaining over MC but at present this is not needed.
    old_sim_df = load_csv(path,old_path=old_dat_path,new_path=old_sim_path,aliases=aliases)

    # Get counts OR y errors from csv
    new_sim_df_count = new_sim_df[count_key] #if yerr_key is None or yerr_key == '' else 1.0/np.square(new_sim_df[yerr_key])
    old_sim_df_count = old_sim_df[count_key] #if yerr_key is None or yerr_key == '' else 1.0/np.square(old_sim_df[yerr_key])
    old_dat_df_count = old_dat_df[count_key] #if yerr_key is None or yerr_key == '' else 1.0/np.square(old_dat_df[yerr_key])

    # Compute scaled quantities
    yerrs            = old_dat_df[yerr_key]
    acceptanceratio  = np.divide(new_sim_df_count,old_sim_df_count) / xs_ratio
    acceptanceratio[np.isinf(acceptanceratio)] = 0
    acceptanceratio[np.isnan(acceptanceratio)] = 0 #NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    scaling          = acceptanceratio * lumi_ratio
    new_dat_df_count = np.multiply(scaling,old_dat_df_count)
    err_scaling      = np.sqrt(np.divide(old_dat_df_count,new_dat_df_count)) #NOTE: SCALE ERRORS ASSUMING POISSONIAN STATISTICS -> d ~ 1/sqrt(N) -> MULTIPLY BY sqrt(N_old_data/N_new_data)
    err_scaling[np.isinf(err_scaling)] = 0
    err_scaling[np.isnan(err_scaling)] = 0 #NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    scaled_yerrs     = np.multiply(err_scaling,yerrs) * 1.0/(tpol_factor * tdil_factor)

    # Set y values to constant and update scaled y errors if requested
    scaled_ys = old_dat_df[y_key] if yvalue<-1 else [yvalue for i in range(len(old_dat_df[y_key]))]
    if yvalue>=-1: scaled_yerrs *= np.sqrt(1-np.square(yvalue*tpol_factor))

    # Copy the old dataframe into the new dataframe
    new_dat_df = old_dat_df.copy(deep=True)

    # Set the entries of your new dataframe
    new_dat_df[count_key] = new_dat_df_count
    new_dat_df[y_key] = scaled_ys
    new_dat_df[yerr_key] = scaled_yerrs
    new_dat_df['scaling'] = scaling
    new_dat_df['acceptanceratio'] = acceptanceratio

    # Save the new dataframe to CSV
    outpath = path.replace('.csv','_rescaled.csv') if len(outpath)==0 else outpath
    new_dat_df.to_csv(outpath, float_format=float_format, index=False)

def get_cut_array(
        cut_titles,
        proj_ids,
        arr_vars,
    ):
    """
    Parameters
    ----------
    cuts : dict, required
        Dictionary of bin ids to list of bin cuts by variable
    proj_ids : list, required
        Array of unique integer bin ids
    arr_vars : list, required
        Variables in which to construct a grid of cuts

    Returns
    -------
    list
        An array of bin cut dictionaries with the same grid structure as :obj:`proj_ids`.

    Raises
    ------
    TypeError
        Raise an error the shape of the array of bin indices is not in :math:`(1,2)`.

    Description
    -----------
    Create a grid of dictionaries of array variables to bin cut titles given a
    dictionary of bin ids to cut title lists and a grid array of projection bin indices.
    Note that the returned array will have the same shape as the given :obj:`proj_ids`.
    """

    # Get grid shape from projection ids array
    shape = np.shape(proj_ids)

    # Create a graph array in the 2D grid case
    if len(shape)==3:
        return [[
            cut_titles[proj_ids[i][j][0]] #NOTE: All BINS IN A 1D BINNING PROJECTION SHOULD HAVE THE SAME BIN CUTS FOR THE ARRAY VARIABLES.
            for j in range(shape[1])] for i in range(shape[0])
        ]

    # Create a graph array in the 1D grid case
    elif len(shape)==2:
        return [
            cut_titles[proj_ids[i][0]] #NOTE: All BINS IN A 1D BINNING PROJECTION SHOULD HAVE THE SAME BIN CUTS FOR THE ARRAY VARIABLES.
            for i in range(shape[0])
        ]

    # Raise an error if another shape length is encountered
    else:
        raise TypeError('`get_cut_array` : `proj_ids` must have len(shape) in (2,3) but shape = ',shape)

def add_cut_array(
        args_array,
        cut_array,
        arr_vars,
    ):
    """
    Parameters
    ----------
    args_array : list, required
        Array of argument dictionaries
    cut_array : list, required
        List array of dictionaries of array variables to bin cut titles
    proj_ids : list, required
        Array of unique integer bin ids
    arr_vars : list, required
        Variables in which to construct a grid of cuts

    Returns
    -------
    list
        An array of bin cut dictionaries with the same grid structure as :obj:`proj_ids`.

    Raises
    ------
    TypeError
        Raise an error the shape of the arrays do not match or the length of the shapes is not in :math:`(1,2)`.

    Description
    -----------
    Add bin cut :obj:`title` and :obj:`ylabel` arguments to each argument dictionary in a grid array
    from a cut array produced from :meth:`get_cut_array`.
    Note that the returned array will have the same shape as the given :obj:`args_array`.
    """

    # Get grid shape from projection ids array
    shape = np.shape(args_array)

    # Check shapes match
    if shape!=np.shape(cut_array): raise TypeError('`add_cut_array` : `args_array` must have same shape as `cut_array` but shapes = ',shape,np.shape(cut_array))

    # Create a graph array in the 2D grid case
    if len(shape)==2: #NOTE: Since this is a 2D grid of dictionaries the numpy shape will only have length 2
            for i in range(shape[0]):
                for j in range(shape[1]):
                    args_array[i][j]['title'] = cut_array[i][j][arr_vars[0]]
                    args_array[i][j]['ylabel'] = cut_array[i][j][arr_vars[1]]

    # Create a graph array in the 1D grid case
    elif len(shape)==1:
        for i in range(shape[0]):
            args_array[i]['title'] = cut_array[i][arr_vars[0]]

    # Raise an error if another shape length is encountered
    else:
        raise TypeError('`add_cut_array` : `args_array` must have len(shape) in (2,3) but shape = ',shape)

def get_bin_mig_mat(
        df,
        id_gen_key='binid_gen',
        id_rec_key='binid_rec',
        mig_key='mig',
    ):
    """
    Parameters
    ----------
    df : pandas.DataFrame, required
        Pandas dataframe of bin migration matrix
    id_gen_key : str, optional
        Key for generated bin indices
    id_rec_key : str, optional
        Key for reconstructed bin indices
    mig_key : str, optional
        Key for bin migration fractions

    Returns
    -------
    np.array
        Bin migration matrix whose indices map generated to reconstructed bins.

    Raises
    ------
    TypeError
        Raise an error if the generated and reconstructed unique bin ids do not have the same members.

    Description
    -----------
    Create a 2D bin migration matrix from a dataframe of generated and reconstructed bin indices
    and the bin migration fractions.
    """

    # Get unique bin ids and check they match between generated and reconstructed
    unique_binids_gen = sorted(df[id_gen_key].unique().tolist())
    unique_binids_rec = sorted(df[id_rec_key].unique().tolist())
    if np.any([el not in unique_binids_gen for el in unique_binids_rec]) or np.any([el not in unique_binids_rec for el in unique_binids_gen]):
        raise TypeError(f"`get_bin_migration_matrix` : Generated and reconstructed unique bin ids must have same members but found: unique_binids_gen={unique_binids_gen}, unique_binids_rec={unique_binids_rec}")

    # Now reshape your bin migration matrix
    shape = [len(unique_binids_gen) for i in range(2)] #NOTE: THIS MUST BE A SQUARE MATRIX
    mig_array = np.reshape(df[mig_key].tolist(),shape)

    return mig_array

def get_subset(
        df,
        bin_ids,
        id_key='bin_id'
    ):
    """
    Parameters
    ----------
    df : pandas.DataFrame, required
        Pandas dataframe containing bin ids and data

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe subset containing only elements whose :obj:`id_key` entries are in :obj:`bin_ids`

    Description
    -----------
    Load a subset of a dataset.
    """
    return df.loc[[el in bin_ids for el in df[id_key]]]

def offset_graph_x(
        g,
        offset,
        axis=0
    ):
    """
    Parameters
    ----------
    g : array, required
        Array describing a 1D graph with structure :obj:`graph[{x_mean,x_err,y_mean,y_err,...},{nbins}]`
    offset : float, required
        Value by which to offset graph values
    axis : int, optional
        Axis along which to offset the graph

    Description
    -----------
    Offset the x values of a graph.
    """

    for idx in range(len(g[axis])):
        g[axis][idx] += offset

def save_txt(
        filename,
        data,
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    ):
    """
    Parameters
    ----------
    filename : str, required
        Output file name
    data : array, required
        2D data array with dimensions :obj:`[N_COLUMNS,N_ROWS]`
    delimiter : str, optional
        CSV format delimiter
    header : str, optional
        CSV header
    fmt : str, optional
        CSV column formats
    comments : str, optional
        CSV comments

    Description
    -----------
    Save a square data array of dimensions :obj:`[N_COLUMNS,N_ROWS]` to a text file.
    """

    # Save to CSV
    if header is None: header = ' '+delimiter+delimiter.join([str(i+1) for i in range(len(data))])#NOTE: ASSUME DATA HAS DIMENSION: [NCOL,NROWS]
    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments)

def save_graph_to_csv(
        filename,
        ct,
        x,
        y,
        xerr=None,
        yerr=None,
        xerr_syst=None,
        yerr_syst=None,
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    ):
    """
    Parameters
    ----------
    filename : str, required
        Output file name
    ct : list, required
        Graph count values with shape :obj:`(nbins)`
    x : list, required
        Graph x values with shape :obj:`(nbins)`
    y : list, required
        Graph y values with shape :obj:`(nbins)`
    xerr : list, optional
        Graph x error values with shape :obj:`(nbins)`
    yerr : list, optional
        Graph y error values with shape :obj:`(nbins)`
    xerr_syst : list, optional
        Graph x systematic error values with shape :obj:`(nbins)` or :obj:`(nbins,2)`
    yerr_syst : list, optional
        Graph y systematic error values with shape :obj:`(nbins)` or :obj:`(nbins,2)`
    delimiter : str, optional
        CSV format delimiter
    header : str, optional
        CSV header
    fmt : str, optional
        CSV column formats
    comments : str, optional
        CSV comments

    Raises
    ------
    ValueError
        Raise an error if the shape of the systematics arrays do not match :obj:`(nbins)` or :obj:`(nbins,2)`.


    Description
    -----------
    Write a graph to a CSV file with optional errors and systematic errors.
    Systematic errors may have high and low values.
    """

    # Create data array
    data = []
    if ct is None or len(ct)==0: ct = [0.0 for el in x]
    if xerr is None or len(xerr)==0: xerr = [0.0 for el in x]
    if yerr is None or len(yerr)==0: yerr = [0.0 for el in x]
    if xerr_syst is None or len(xerr_syst)==0: xerr_syst = [0.0 for el in x]
    if yerr_syst is None or len(yerr_syst)==0: yerr_syst = [0.0 for el in x]
    xerr_syst_shape = np.shape(xerr_syst)
    yerr_syst_shape = np.shape(yerr_syst)
    for i, el in enumerate(x):
        data_i = [i, ct[i], x[i], y[i], xerr[i], yerr[i]]
        if len(xerr_syst_shape)==1 and len(yerr_syst_shape)==1:
            data_i.extend([xerr_syst[i], yerr_syst[i]])
        elif len(xerr_syst_shape)==2 and xerr_syst_shape[1]==2 and len(yerr_syst_shape)==1:
            data_i.extend([xerr_syst[i][0], xerr_syst[i][1], yerr_syst[i]])
        elif len(xerr_syst_shape)==1 and len(yerr_syst_shape)==2 and yerr_syst_shape[1]==2:
            data_i.extend([xerr_syst[i], yerr_syst[i][0], yerr_syst[i][1]])
        elif len(xerr_syst_shape)==2 and xerr_syst_shape[1]==2 and len(yerr_syst_shape)==2 and yerr_syst_shape[1]==2:
            data_i.extend([xerr_syst[i][0], xerr_syst[i][1], yerr_syst[i][0], yerr_syst[i][1]])
        else:
            raise ValueError(f"ERROR: xerr_syst_shape={xerr_syst_shape} or yerr_syst_shape={yerr_syst_shape} does not have shape (nbins) or (nbins,2).")
        data.append(data_i)
    data = np.array(data)

    # Save data to file
    save_txt(filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments)

def save_graph_systematics_to_csv(
        filename,
        x,
        yerrs_syst=None,
        delimiter=",",
        header=None,
        fmt=None,
        comments=''
    ):
    """
    Parameters
    ----------
    filename : str, required
        Output file name
    x : list, required
        Graph x values with shape :obj:`(nbins)`
    yerr_syst : list, optional
        Graph y systematic error values decomposed into the different sources of systematic error with shape :obj:`(nbins,nsources)` or :obj:`(nbins,nsources,2)`
    delimiter : str, optional
        CSV format delimiter
    header : str, optional
        CSV header
    fmt : str, optional
        CSV column formats
    comments : str, optional
        CSV comments

    Raises
    ------
    ValueError
        Raise an error if the shape of the systematics is does not match :obj:`(nbins,nsources)` or :obj:`(nbins, nsources,2)`.

    Description
    -----------
    Write a set of graph y systematic errors to a CSV file with the systematic error values broken down by source and allowing high and low errors.
    This means the argument :obj:`yerr_syst` should have shape :obj:`(nbins, nsources)` or :obj:`(nbins, nsources,2)`
    where :obj:`nbins` is the number of kinematic bins and :obj:`nsources` is the number of sources of systematic error.
    """

    # Create data array
    data = []
    if yerrs_syst is None or len(yerrs_syst)==0: yerrs_syst = [[0.0] for el in x]
    yerrs_syst_shape = np.shape(yerrs_syst)        
    if yerrs_syst_shape[0]==len(x) and len(yerrs_syst_shape)==2:
        for i, el in enumerate(x):
            data.append([i, x[i], *yerrs_syst[i]])
    elif yerrs_syst_shape[0]==len(x) and (len(yerrs_syst_shape)==3 and yerrs_syst_shape[2]==2):
        for i, el in enumerate(x):
            data_i = [i, x[i]]
            for source in yerrs_syst[i]:
                data_i.extend(*source)
            data.append(data_i)
    else:
        raise ValueError(f"ERROR: yerrs_syst has shape {yerrs_syst} but allowed shapes are ({len(x)},*) and ({len(x)},*,2).")
    data = np.array(data)

    # Save data to file
    save_txt(filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments)

def save_bin_mig_mat_to_csv(
        bin_mig_mat,
        base_dir='./',
        basename='',
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    ):
    """
    Parameters
    ----------
    bin_mig_mat : np.array, required
        2D bin migration matrix with :obj:`(i,j) -> (generated,reconstructed)`
    base_dir : str, required
        Path to directory in which matrix will be saved
    basename : str, optional
        Name of reconstructed bin variable
    delimiter : str, optional
        CSV format delimiter
    header : str, optional
        CSV header
    fmt : str, optional
        CSV column formats, automatically set if not specified
    comments : str, optional
        CSV comments

    Raises
    ------
    TypeError
        Raise an error if the bin migration matrix is not square

    Description
    -----------
    Save a 2D bin migration matrix mapping generated bins to reconstructed bins to a CSV file
    with an added initial row and column for the bin indices.  Note that files will be saved
    to :obj:`<base_dir>/bin_mig_mat_<basename>.csv`.
    """

    # Check bin migration matrix shape
    if np.shape(bin_mig_mat)[0]!=np.shape(bin_mig_mat)[1] or len(np.shape(bin_mig_mat))!=2:
        raise TypeError("Bin migration matrix must be square but has shape "+str(np.shape(bin_mig_mat)))

    # Set output filename
    filename = 'bin_mig_mat_'+basename+'.csv'
    filename = os.path.join(base_dir,filename)

    # Create new table with int bin labels
    nbins = np.shape(bin_mig_mat)[0]
    new_shape = list(np.shape(bin_mig_mat)) #NOTE: List is important here!
    new_shape[1] += 1 #NOTE: Increase the number of columns to accomodate bin numbers in the initial column
    data = np.zeros(new_shape)
    data[:,0] = [i for i in range(1,nbins+1)]
    data[0:,1:] = bin_mig_mat

    # Set column formats if not given
    if fmt is None:
        fmt = ["%.3g" for i in range(np.shape(bin_mig_mat)[0])]
        fmt = ["%d",*fmt]

    save_txt(filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments)

def apply_bin_mig(
        df,
        inv_bin_mig_mat,
        results_keys = ['a0'],
    ):
    """
    Parameters
    ----------
    df : pandas.DataFrame, required
        Pandas dataframe of asymmetry results
    inv_bin_mig_mat : np.array, required
        Inverse of bin migration matrix mapping generated bins to reconstructed bins
    results_keys : list, optional
        List of keys to results entries to which to apply bin migration corrections

    Description
    -----------
    Apply a bin migration correction to a set of results contained within a pandas dataframe.
    """

    # Then multiply results and inverse bin migration and reset the original dataframe
    for result_key in results_keys:
        df[result_key] = np.matmul(inv_bin_mig_mat,df[result_key])

def compute_systematics(
        results,
        bin_migration_mat=None,
        systematic_scales_mat=None,
        systematics_additive_mat=None
    ):
    """
    Parameters
    ----------
    results : np.array, required
        Array of asymmetry results
    bin_migration_mat : np.array, optional
        Bin migration matrix mapping generated bins to reconstructed bins
    systematic_scales_mat : np.array, optional
        Array of systematic errors to be scaled by y values before being added in quadrature to other systematic errors
    systematics_additive_mat : np.array, optional
        Array of absolute systematic errors to add to other systematic errors

    Returns
    -------
    np.array
        Array of systematics values added in quadrature

    Description
    -----------
    Compute the systematic errors for a 1D binning scheme given any combination
    of bin migration matrix (this will be inverted internally with :meth:`np.linalg.inv`),
    a systematic scaling array to multiply by the actual results and be added in quadrature,
    and a systematic additive array to add to the other systematic errors **not** in quadrature.
    """

    # Initialize empty systematics array
    systematics = np.zeros(results.shape)

    # Compute and add bin migration systematics
    if bin_migration_mat is not None: #NOTE: ASSUME BINNING IS 1D HERE.
        bin_migration_mat_inv = np.linalg.inv(bin_migration_mat) #NOTE: Make sure that bin_migration_mat is of the form f[i,j] = [# gen in bin i AND rec. in bin j] / [# gen. in bin i], remember that ROOT histogram indexing i,j is opposite (idx_x,idx_y) typical matrix indexing (idx_y,idx_x) and begins at 1 not 0
        new_systematics = np.add(results,-np.matmul(bin_migration_mat_inv,results)) # DeltaA = a - f_inv . a
        systematics = np.sqrt(np.square(systematics) + np.square(new_systematics))

    # Apply multiplicative scale systematics, note that these should already be summed over all sources of systematic error
    if systematic_scales_mat is not None:
        systematics = np.sqrt(np.square(systematics) + np.square(np.multiply(results,systematic_scales_mat))) #NOTE: IMPORTANT!  ADD IN QUADRATURE.

     # Apply additive scale systematics, note that these should already be summed over all sources of systematic error
    if systematics_additive_mat is not None:
        systematics += systematics_additive_mat

    return systematics

def set_default_plt_settings():
    """
    Description
    -----------
    Set plt.rc parameters for font sizes and family and tick font size and tick length and direction
    in a nice format.
    """

    # Set font sizes
    plt.rc('font', size=25) #controls default text size
    plt.rc('axes', titlesize=50) #fontsize of the title
    plt.rc('axes', labelsize=50) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=25) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=25) #fontsize of the y tick labels
    plt.rc('legend', fontsize=20) #fontsize of the legend

    # Get some nicer plot settings
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['figure.autolayout'] = True

    # Set tick parameters
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)

def plot_injected_asyms(
        ax1,
        asyms,
        ytitles,
        colors,
        sgasym_idx = 0,
        ylims = (-1.0,1.0),
        label_base='Injected Signal ',
        linestyle='--',
        linewidth=1,
    ):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    asyms : list, required
        List of injected asymmetries 
    ytitles : list, required
        List of y axis titles
    colors : list, required
        List of colors for each injected asymmetries
    sgasym_idx : int, optional
        Injected signal asymmetry index
    ylims : tuple, optional
        y limits for plotting
    label_base : str, optional
        Base label to prepend to ytitles for injected asymmetries
    linestyle : str, optional
        Line style
    linewidth : int, optional
        Line width

    Description
    -----------
    Plot the injected asymmetries for each bin in a 1D binning scheme offsetting repeat values
    by a small amount.  Note that injected asymmetries should have shape :obj:`(N_ASYMS)` if plotting
    constant asymmetries or :obj:`(N_ASYMS,2,N_POINTS)` if you would like to plot function data :obj:`(x,y)`.
    """

    # Loop asymmetries and plot using a small offset for constant asymmetries
    plotted_values = {} #NOTE: Keep track of how many times you've plotted each constant asymmetry value
    for idx in range(len(asyms)):

        # Plot injected asymmetries as (x,y) data OR axis lines
        if len(np.shape(asyms))>1:
            ax1.plot(asyms[idx][0], asyms[idx][1], color=colors[idx] if idx<len(colors) else None, linestyle=linestyle, linewidth=linewidth, alpha=0.5 if idx!=sgasym_idx else 1.0, label=label_base+ytitles[idx] if idx<len(ytitles) else None)
        else:
            # Set the offset
            offset = 0.0
            if asyms[idx] in plotted_values.keys():
                offset = plotted_values[asyms[idx]] * 0.0025 * (ylims[1] - ylims[0])
                plotted_values[asyms[idx]] += 1
            else:
                plotted_values[asyms[idx]] = 1
            
            # Flip offset if the asymmetry is negative
            if asyms[idx]<0.0: offset *= -1.0

            # Plot the constant asymmetry
            ax1.axhline(asyms[idx]+offset, color=colors[idx] if idx<len(colors) else None, linestyle=linestyle, linewidth=linewidth, alpha=0.5 if idx!=sgasym_idx else 1.0, label=label_base+ytitles[idx] if idx<len(ytitles) else None)

def plot_watermark(
        ax1,
        watermark='CLAS12 Preliminary'
    ):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    watermark : str, optional
        Watermark text

    Description
    -----------
    Plot a watermark.
    """
    plt.text(0.5, 0.5, watermark,
                size=50, rotation=25., color='gray', alpha=0.25,
                horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)

def plot_vlines(
        hist,
        binlims = [],
        linestyle = 'dotted',
    ):
    """
    Parameters
    ----------
    hist : tuple, required
        Matplotlib.pyplot histogram of y and x values (np.ndarray, np.ndarray, ...)
    binlims : list, required
        List of bin limits in a 1D binning scheme
    linestyle : str, optional
        Line style

    Description
    -----------
    Draw vertical bin limit lines on a histogram
    """

    # Loop middle bin limits
    for xval in binlims[1:-1]:
        
        # Loop histogram x values
        for idx in range(len(hist[1])-1):

            # Check if bin limit is in bin
            binx_low, binx_high = hist[1][idx], hist[1][idx+1]
            if xval>= binx_low and xval<binx_high:

                # Plot lower bin limit
                ymax = hist[0][idx]
                plt.vlines(xval, 0.0, ymax, linestyle=linestyle)

def plot_hists(
        ax1,
        hist_paths = [],
        hist_keys = [],
        clone_axis = True,
        ylabel = 'Density',
        ylims = (0.0,0.05),
        histtype = 'step',
        hist_colors = None,
        alpha=0.5,
        linewidth=2,
        density=True,
        log=False,
        hist_labels = None,
        binlims = [],
        vlinestyle = 'dotted',
        vline_hist_idx = -1,
        legend_loc = 'upper right',
        hist_dim = 1,
    ):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    clone_axis : bool, optional
        Option to create a twin y axis sharing the x axis of the given axis
    ylabel : str, optional
        Y axis label for twin axis
    ylims : tuple, optional
        Y axis limits for twin axis
    histtype : str, optional
        Matplotlib.pyplot histogram type
    hist_colors : list, optional
        List of histogram colors
    alpha : float, optional
        Alpha plotting paramater for histograms
    linewidth : int, optional
        Line width for plotting histograms
    density : bool, optional
        Option to normalize histograms
    log : bool, optional
        Option to plot y-axis on a log scale
    hist_labels : list, optional
        List of histogram labels
    binlims : list, optional
        List of bin limits in a 1D binning scheme
    vlinestyle : str, optional
        Vertical line style for drawing bin limits on histogram
    vline_hist_idx : int, optional
        Index of histogram for which to draw vertical lines for bin limits
    legend_loc : str, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to :obj:`None` or :obj:`''`
    hist_dim : int, optional
        Dimension of histogram to plot.  Must be either 1 or 2.

    Description
    -----------
    Draw a set of histograms on a matplotlib.pyplot axis, optionally cloning the given axis
    and optionally drawing vertical bin limits on a single histogram.
    """

    # Clone y axis and set labels
    ax2 = ax1.twinx() if clone_axis else ax1  # instantiate a second y-axis that shares the same x-axis
    if clone_axis:
        ax2.set_ylabel(ylabel)
        ax2.set_ylim(*ylims)

    # Set colors to just be taken from the current color palette if not provided
    if hist_colors is None:
        hist_colors = [None for i in range(len(hist_paths))]

    # Set histogram labels to just be taken from the current color palette if not provided
    if hist_labels is None:
        hist_labels = ["h"+str(i) for i in range(len(hist_paths))]

    # Check line width as a flag for drawing th2ds
    if hist_dim==1:

        # Loop histograms and draw
        for idx, hist_path in enumerate(hist_paths):

            # Load histogram and convert to numpy
            h_y, h_bins = load_th1(hist_path,hist_keys[idx])

            # Get mean x bin values
            h_x = [(h_bins[i]+h_bins[i+1])/2 for i in range(len(h_bins)-1)]

            # Plot histogram
            h = ax2.hist(
                h_x,
                bins=h_bins,
                weights=h_y/np.sum(h_y) if density else h_y,
                histtype=histtype,
                color=hist_colors[idx],
                alpha=alpha,
                linewidth=linewidth,
                label=hist_labels[idx],
                density=False,
                log=log,
            )

            # Plot bin limits if supplied and we are on the last histogram
            if idx==(vline_hist_idx if vline_hist_idx>=0 else len(hist_paths)-1) and len(binlims)>0:
                plot_vlines(
                    h,
                    binlims,
                    linestyle = vlinestyle,
                )

    # Assume TH2D histograms and plot
    else:

        # Loop histograms and draw
        for idx, hist_path in enumerate(hist_paths):

            # Load histogram and convert to numpy
            h2 = load_th1(hist_path,hist_keys[idx])

            # Plot histogram
            plot_th2(h2, ax1, add_colorbar=True, norm=colors.LogNorm(), label=hist_labels[idx]) #NOTE: Just fix the colorbar option for now since it's not likely you would want a 2D hist without a z axis scale legend.

    # Plot legend if you cloned axis
    if clone_axis and legend_loc is not None and legend_loc!='': ax2.legend(loc=legend_loc)

def get_bin_kinematics_title(
        bin_id,
        df,
        cols=['x','Q2'],
        col_titles={'x':'x', 'Q2':'Q^{2}'},
        err_ext='_err',
        sep=' , '
    ):
    """
    Parameters
    ----------
    bin_id : int, required
        Bin id in kinematics dataframe
    df : pandas.Dataframe, required
        Dataframe of kinematic variable means and errors in each kinematic bin
    cols : list, optional
        List of column names of kinematics to add to bin title
    col_titles : list, optional
        List of kinematics LaTeX titles
    err_ext : str, optional
        Extension for forming column names of kinematic variable errors
    sep : str, optional
        Separator string for kinematics values in bin title

    Returns
    -------
    str
        Title string

    Description
    -----------
    Create a title string for the requested bin index and kinematic variables showing the mean and error values for
    each of the requested kinematic variables in the given bin.
    """

    return sep.join(
        [
            f"$<{col_titles[col]}> = {df.iloc[bin_id].loc[col]:.2f}\\pm{df.iloc[bin_id].loc[cols[idx]+err_ext]:.2f}$"
            for idx, col in enumerate(cols)
        ]
    )

def get_lims_coords(
        node,
        outer_xlims,
        outer_ylims,
        var_keys=[],
        nested_key='nested',
        lims_key='lims',
        swap_axes=False,
    ):
    """
    Parameters
    ----------
    node : dict, required
        Bin scheme node
    outer_xlims : list, required
        List of outer limits for the x-axis variable
    outer_ylims : list, required
        List of outer limits for the y-axis variable
    var_keys : list, optional
        List names of x and y axis variables for a grid bin scheme
    nested_key : str, optional
        Key for nested bins
    lims_key : str, optional
        Key for bin limits
    swap_axes : bool, optional
        Swap the default x and y axes order

    Raises
    ------
    ValueError
        Raises error if no acceptable bin scheme is found.

    Returns
    -------
    list
        List of line coordinates in the form :math:`((x_1,x_2),(y_1,y_2))`

    Description
    -----------
    Get a list of line coordinates delineating the bin limits for a nested 2D binning scheme.
    """

    # Initialize coordinates list
    lims_coords = []

    # Check node
    if (type(node)==dict and nested_key in node and type(node[nested_key])==list and type(node[nested_key][0])==dict):

        # Get first level nested node and get horizontal limits coordinates
        binvar_x = list(node[nested_key][0].keys())[0]
        node_nested = node[nested_key][0][binvar_x]
        horizontal_lims = [
            [outer_xlims,[y0,y0]] for y0 in node_nested[lims_key][1:-1]
        ]
        ylims = node_nested[lims_key]

        # Check nested node and get limits list
        xlims = []
        if (type(node_nested)==dict and nested_key in node_nested and type(node_nested[nested_key])==list and type(node_nested[nested_key][0])==dict):

            # Loop nested nodes and get limits lists
            for el in node_nested[nested_key]:
                binvar_y = list(el.keys())[0]
                if lims_key in el[binvar_y]: xlims.append(el[binvar_y][lims_key])

            # Loop xlims and set vertical limit coordinates
            vertical_lims = []
            for yidx, xlim in enumerate(xlims):
                for xidx, x in enumerate(xlim[1:-1]):
                    el = [[xlims[yidx][xidx+1],xlims[yidx][xidx+1]],[ylims[yidx] if yidx>0 else outer_ylims[0],ylims[yidx+1] if yidx<len(xlims)-1 else outer_ylims[-1]]]
                    vertical_lims.append(el)


            # Add to coordinates list
            lims_coords.extend(vertical_lims)
        lims_coords.extend(horizontal_lims)
    
    # Grid scheme case
    elif type(node)==dict and len(var_keys)==2:

        # Get limits from node
        xvar, yvar = var_keys
        xlims = node[xvar]
        ylims = node[yvar]
        
        # Get line coordinates
        horizontal_lims = [
            [outer_xlims,[y0,y0]] for y0 in ylims[1:-1]
        ]
        vertical_lims = [
            [[x0,x0],outer_ylims] for x0 in xlims[1:-1]
        ]

        # Add to coordinates list
        lims_coords.extend(vertical_lims)
        lims_coords.extend(horizontal_lims)

    # Default to raising value error
    else:
        raise ValueError(f'Could not identify a 2D nested or grid bin scheme in `node`:\n{node}')

    # Swap axes if needed
    if swap_axes:
        lims_coords = [[el[1],el[0]] for el in lims_coords] 

    return lims_coords

def plot_lines(
        ax,
        coordinates,
        linecolor='red',
        linewidth=1
    ):
    """
    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot axis to plot on
    coordinates : list of lists, required
        List of coordinate pairs each with structure :math:`((x_1,x_2),(y_1,y_2))`
    color : str, optional
        Line color
    linewidth : int, optional
        Line width

    Description
    -----------
    Plot a set of lines from a list of coordinates.
    """

    # Check coordinates shape
    if len(np.shape(coordinates))==3 and np.shape(coordinates)[1:]!=(2,2):
        raise ValueError(f'Expected shape (2,2) but got {np.shape(coordinates)}')

    # Plot lines
    for coords in coordinates:
        ax.plot(*coords, color=linecolor, linewidth=linewidth, marker = 'o', markersize=0)

def get_bin_centers(
        cuts,
        swap_axes=False
    ):
    """
    Parameters
    ----------
    cuts : dict, required
        Dictionary of bin ids to bin cuts
    swap_axes : bool, optional
        Swap the default x and y axes order

    Returns
    -------
    tuple
        Tuple of maps of bin ids to bin centers and bin widths respectively

    Description
    -----------
    Get the bin center coordinates from a list of (rectangular 2D) bin cuts.
    Note that bin cuts are assumed to be of the form
    :obj:`(binvar_x>=x_min && binvar_x<=x_max) && (binvar_y>=y_min && binvar_y<=y_max)`.
    """

    # Convert dictionary to lists
    new_cuts = [cuts[idx] for idx in cuts]
    bin_ids = [idx for idx in cuts]

    # Convert 2D bin cuts into bin widths and centers
    bin_centers = [cut.replace('=','').split(") && (") for cut in new_cuts]
    bin_centers = [[el[0].replace('(','').split(' && '),el[1].replace(')','').split(' && ')] for el in bin_centers]
    bin_centers = [[[float(el[idx][0].split('>')[1]),float(el[idx][1].split('<')[1])] for idx, _ in enumerate(el)] for el in bin_centers]
    bin_widths  = [[el[0][1]-el[0][0],el[1][1]-el[1][0]] for el in bin_centers] #NOTE: ORDERING IS IMPORTANT HERE
    bin_centers = [[np.average(el[0]).item(),np.average(el[1]).item()] for el in bin_centers]

    # Swap axes if needed
    if swap_axes:
        bin_centers = [[el[1],el[0]] for el in bin_centers]
        bin_widths = [[el[1],el[0]] for el in bin_widths]

    # Convert back to dictionaries
    bin_centers = {bin_ids[idx]:bin_centers[idx] for idx in range(len(bin_centers))}
    bin_widths  = {bin_ids[idx]:bin_widths[idx]  for idx in range(len(bin_widths))}

    return bin_centers, bin_widths

def plot_bin_ids(
        ax,
        bin_centers,
        bin_widths=None,
        size=25,
        color='red',
        alpha=1.0,
    ):
    """
    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot axis to plot on
    bin_centers : dict, required
        Dictionary of bin ids to bin centers
    bin_widths : dict, optional
        Dictionary of bin ids to bin widths each of the form :obj:`(width_x, width_y)`
    size : int, optional
        Font size for bin id text
    color : str, optional
        Text color
    alpha : float, optional
        Text alpha value
    """

    # Plot bin ids on bin centers
    for bin_id in bin_centers:
        ax.text(*bin_centers[bin_id], f'{bin_id}', size=size, color=color, alpha=alpha, horizontalalignment='center', verticalalignment='center')

def plot_th2(
        h2,
        ax,
        add_colorbar=True,
        norm=colors.LogNorm(),
        **kwargs
    ):
    """
    Parameters
    ----------
    h2 : tuple or list, required
        List of 2D histogram data with structure :obj:`(weights, xbins, ybins)`
    ax : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot axis to plot on
    add_colorbar : bool, optional
        Add a colorbar to show the z-axis scale
    norm : str or matplotlib.colors.Normalize, optional
        Normalization used to scale data to :math:`[0,1]` range before mapping to a color map
    **kwargs
        Additional parameters are passed along to :meth:`matplotlib.pyplot.hist2d`

    Description
    -----------
    Easily plot a :obj:`TH2` histogram loaded from ROOT.
    """

    # Get the middle values of each bin
    x = np.ravel([[np.average([h2[1][i],h2[1][i+1]]) for j in range(len(h2[2])-1) ] for i in range(len(h2[1])-1)])
    y = np.ravel([[np.average([h2[2][j],h2[2][j+1]]) for j in range(len(h2[2])-1) ] for i in range(len(h2[1])-1)])

    # Get the counts in each bin
    weights = np.ravel(h2[0])

    # Get the bin sizes
    bins = (h2[1], h2[2])

    # Plot the histogram
    hist2d = ax.hist2d(x,y,weights=weights,bins=bins, norm=norm, **kwargs)
    if add_colorbar: plt.colorbar(hist2d[3],ax=ax)

def plot_systematics(
        x_means,
        yerr_syst,
        palette = 'Dark2',
        stacked = False,
        syst_names = None,
        syst_labels = None,
        xlims   = (0.0,1.0),
        ylims   = (-1.0,1.0),
        xvar    = 'x',
        title   = 'Systematic Errors',
        xtitle  = '$Q^{2} (GeV^{2})$',
        ytitle  = '$\\Delta \\mathcal{A}$',
        outpath = 'systematics.pdf',
        watermark = 'CLAS12 Preliminary',
        use_default_plt_settings = True,
        legend_loc = 'upper left',
        ecolor = 'black',
        elinewidth = 2.0,
        capsize = 18,
        capthick = 2.0,
        marker = 'o',
        markersize = 20,
        linestyle = None,
        linewidth = 0.0,
        gridlinewidth = 0.5,
        axlinewidth = 1.0,
        log = False,
        figsize = (16,10),
    ):
    """
    Parameters
    ----------
    x_means : list, required
        Mean x values for each bin
    yerr_syst : np.array, required
        Array of absolute systematic error in each bin further indexed by the sources of systematic error
    palette : str, optional
        Seaborn color palette
    stacked : bool, optional
        Whether to stack histograms from different sources of systematic error
    syst_names : list, optional
        List of column names for each source of systematic error
    syst_labels : list, optional
        List of labels for each source of systematic error
    xlims : tuple, optional
        x limits for plotting
    ylims : tuple, optional
        y limits for plotting
    xvar : str, optional
        Bin variable name
    title : str, optional
        Plot title
    xtitle : str, optional
        x axis title
    ytitle : str, optional
        y axis title
    outpath : str, optional
        Name of output pdf
    watermark : str, optional
        Optional watermark to put on top of plot
    use_default_plt_settings : bool, optional
        Option to use default font and tick parameter style settings
    legend_loc : str, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to :obj:`None` or :obj:`''`
    ecolor : str, optional
        Error line color
    ecolor : float, optional
        Error line width
    capsize : int, optional
        Error cap size
    capthick : float, optional
        Error cap thickness
    marker : str, optional
        Marker type
    markersize : int, optional
        Marker size
    linestyle : str, optional
        Line style
    linewidth : float, optional
        Line width
    gridlinewidth : float, optional
        Grid line width
    axlinewidth : float, optional
        Axis line and injected asymmetries line width
    log : bool, optional
        Option to plot y-axis on a log scale
    figsize : tuple, optional
        Figure size

    Description
    -----------
    Plot the systematic error for each bin in a 1D binning scheme broken down by sources of systematic error.
    Save systematics breakdowns to CSV in :obj:`<outpath>.csv`.  Note that this does **not** allow for asymmetric errors.
    """

    # Set color palette
    sbn.set_palette(palette)

    # Use default plotting settings
    if use_default_plt_settings: set_default_plt_settings()

    # Set up plot
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)

    # Plot systematics by source for each x point
    nbins = len(x_means)
    xbins = np.moveaxis(np.array([x_means for el in range(np.shape(yerr_syst)[1])]),(0,1),(1,0))
    s1 = plt.hist(xbins, weights=yerr_syst, bins=nbins, alpha=0.5, label=syst_labels, stacked=stacked, log=log)

    # Plot zero line
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)

    # Add water mark
    if watermark is not None and watermark!='': plot_watermark(ax1,watermark=watermark)

    # Plot legend
    if legend_loc is not None and legend_loc!='': ax1.legend(loc=legend_loc)

    # Save figure
    f1.savefig(outpath)

    # Save plot data to csv
    delimiter = ","
    if syst_names is None: syst_names = ["syst"+str(idx) for idx in range(nbins)]
    header    = delimiter.join(["bin",*syst_names]) #NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    syst_fmts = ["%.3g" for idx in range(nbins)]
    fmt       = ["%d",*syst_fmts]
    comments  = ""

    # Save to CSV
    save_graph_systematics_to_csv(
        outpath+'.csv',
        x_means,
        yerrs_syst=yerr_syst,
        delimiter=delimiter,
        header=header,
        fmt=fmt,
        comments=comments
    )

def plot_results(
        ax1,
        ct_mean = None,
        x_mean = None,
        y_mean = None,
        xerr_mean = None,
        yerr_mean = None,
        xerr_syst = None,
        yerr_syst = None,
        y_min  = None,
        y_max  = None,
        y_std  = None,
        ydiff_mean = None,
        ydiff_std = None,
        ydiff_min = None,
        ydiff_max = None,
        xlims = [0.0,1.0],
        ylims = [0.0,1.0],
        title = 'Asymmetry Results',
        xvar  = 'x',
        xlabel = '$x$',
        ylabel = '$\\mathcal{A}$',
        sgasym_labels = ['$\\mathcal{A}$'],
        bgasym_labels = ['$\\mathcal{A}$'],
        sgasym_idx = 0,
        sgasyms = [0.10],
        bgasyms = [0.00],
        sg_colors  = ['blue'],
        bg_colors  = ['red'],
        fill_color = 'gray',
        outpath = 'out.pdf',
        watermark = 'CLAS12 Preliminary',
        show_injected_asymmetries = False,
        legend_loc = 'upper left',
        ecolor = 'black',
        elinewidth = 2.0,
        capsize = 18,
        capthick = 2.0,
        marker = 'o',
        markersize = 20,
        linestyle = None,
        linewidth = 0.0,
        gridlinewidth = 0.5,
        axlinewidth = 1.0,
        hist_paths = [],
        hist_keys = [],
        hist_clone_axis = True,
        hist_ylabel = 'Density',
        hist_ylims = (0.0,0.05),
        histtype = 'step',
        hist_colors = None,
        hist_alpha=0.5,
        hist_linewidth=2,
        hist_density=True,
        hist_log=False,
        hist_labels = None,
        binlims = [],
        vlinestyle = 'dotted',
        vline_hist_idx = -1,
        hist_legend_loc = 'upper right',
        hist_dim = 1,
        old_dat_path = 'old_dat_path.csv',
        new_sim_path = 'new_sim_path.csv',
        old_sim_path = 'old_sim_path.csv',
        count_key = 'count',
        yerr_key = '',
        xs_ratio = 1.0,
        lumi_ratio = 0.0,
        tpol_factor = 1.0,
        tdil_factor = 1.0,
        graph_yvalue = -100.0,
        aliases = {},
        plot_xerrors = False,
    ):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    ct_mean : list, optional
        Count mean values for each bin with shape :obj:`(nbins)`
    x_mean : list, optional
        x mean values for each bin with shape :obj:`(nbins)`
    y_mean : list, optional
        y mean values for each bin with shape :obj:`(nbins)`
    xerr_mean : list, optional
        x error alues for each bin with shape :obj:`(nbins)`
    yerr_mean : list, optional
        y error values for each bin with shape :obj:`(nbins)`
    xerr_syst : list, optional
        x systematic error alues for each bin with shape :obj:`(nbins)` or :obj:`(nbins,2)`
    yerr_syst : list, optional
        y systematic error values for each bin, with shape :obj:`(nbins)` or :obj:`(nbins,2)`
    y_min : list, optional
        y minimum values for each bin with shape :obj:`(nbins)`
    y_max : list, optional
        y maximum values for each bin with shape :obj:`(nbins)`
    y_std : list, optional
        y standard deviation values for each bin with shape :obj:`(nbins)`
    ydiff_mean : list, optional
        y difference from injected signal asymmetry mean values for each bin with shape :obj:`(nbins)`
    ydiff_std : list, optional
        y difference from injected signal asymmetry standard deviation values for each bin with shape :obj:`(nbins)`
    ydiff_min : list, optional
        y difference from injected signal asymmetry minimum values for each bin with shape :obj:`(nbins)`
    ydiff_max : list, optional
        y difference from injected signal asymmetry maximum values for each bin with shape :obj:`(nbins)`
    xlims : tuple, optional
        x limits for plotting
    ylims : tuple, optional
        y limits for plotting
    title : str, optional
        Plot title
    xvar : str, optional
        Bin variable name
    xlabel : str, optional
        x axis label
    ylabel : str, optional
        y axis label
    sgasym_labels : list, optional
        List of signal asymmetry labels
    bgasym_labels : list, optional
        List of background asymmetry labels
    sgasym_idx : int, optional
        Index of injected signal asymmetry
    sgasyms : list, optional
        List of injected signal asymmetries
    bgasyms : list, optional
        List of injected background asymmetries
    sg_colors : list, optional
        List of signal asymmetry plotting colors
    bg_colors : list, optional
        List of background asymmetry plotting colors
    fill_color : str, optional
        Color of 1 sigma band or systematic uncertainties
    outpath : str, optional
        Name of output pdf
    watermark : str, optional
        Optional watermark to put on top of plot
    show_injected_asymmetries : bool, optional
        Option to show injected signal and background asymmetries
    legend_loc : str, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to :obj:`None` or :obj:`''`
    ecolor : str, optional
        Error line color
    ecolor : float, optional
        Error line width
    capsize : int, optional
        Error cap size
    capthick : float, optional
        Error cap thickness
    marker : str, optional
        Marker type
    markersize : int, optional
        Marker size
    linestyle : str, optional
        Line style
    linewidth : float, optional
        Line width
    gridlinewidth : float, optional
        Grid line width
    axlinewidth : float, optional
        Axis line and injected asymmetries line width
    hist_clone_axis : bool, optional
        Option to create a twin y axis sharing the x axis of the given axis
    hist_ylabel : str, optional
        Y axis label for twin axis
    hist_ylims : tuple, optional
        Y axis limits for twin axis
    histtype : str, optional
        Matplotlib.pyplot histogram type
    hist_colors : list, optional
        List of histogram colors
    hist_alpha : float, optional
        Alpha plotting paramater for histograms
    hist_linewidth : int, optional
        Line width for plotting histograms
    hist_density : bool, optional
        Option to normalize histograms
    hist_log : bool, optional
        Option to plot y-axis on a log scale
    hist_labels : list, optional
        List of histogram labels
    binlims : list, optional
        List of bin limits in a 1D binning scheme
    vlinestyle : str, optional
        Vertical line style for drawing bin limits on histogram
    vline_hist_idx : int, optional
        Index of histogram for which to draw vertical lines for bin limits
    hist_legend_loc : str, optional
        Matplotlib.pyplot legend location string for histograms, will not be plotted if set to :obj:`None` or :obj:`''`
    hist_dim : int, optional
        Dimension of histogram to plot.  Must be either 1 or 2.
    old_dat_path : str, optional
        Part of :obj:`path` to replace
    new_sim_path : str, optional
        Path with which to replace :obj:`old_dat_path` to find new simulation graph CSVs
    old_sim_path : str, optional
        Path with which to replace :obj:`old_dat_path` to find new simulation graph CSVs
    count_key : str, optional
        CSV column key for graph bin counts
    yerr_key : str, optional
        CSV column key for graph bin y errors
    xs_ratio : float, optional
        Cross-section ratio (new/old) for scaling
    lumi_ratio : float, optional
        Luminosity ratio (new/old) for scaling
    tpol_factor : float, optional
        Target polarization factor for rescaling
    tdil_factor : float, optional
        Target dilution factor for rescaling
    graph_yvalue : float, optional
        Constant asymmetry value to be plotted for showing rescaled errors
    aliases : dict, optional
        Map of configuration option names to maps of option values to string aliases for chained CSVs when rescaling
    plot_xerrors : bool, optional
        Option to plot x errors on asymmetry graph

    Description
    -----------
    Plot asymmetry results for each bin in a 1D binning scheme showing projection variable histograms, bin limits, systematic errors,
    a standard deviation band for aggregate graphs, injected asymmetries, watermark, and legend if desired.  Save results
    and differences from the injected signal to CSV in :obj:`<outpath>.csv` and :obj:`<outpath>_ydiff.csv`.
    Optionally, rescale the input graph using :meth:`rescale_graph_data`.  Note that the graph y values can be
    set to a constant with :obj:`graph_yvalue` if you only wish to show the rescaled errors, and that here, the rescaled asymmetry
    errors will be further scaled to account for target polarization and dilution factors like
    :math:`\\sigma_{A} = \\frac{1}{P_{Target} \\cdot D_{Target}}\\cdot\\sqrt{\\frac{1-(A \\cdot P_{Target})^{2}}{N_{Rescaled}}}`.
    """

    # Rescale graph
    scaling, acceptanceratio = None, None
    rescale = lumi_ratio>0.0
    if rescale:

        # Get rescaled graph data
        rescaled_graph = rescale_graph_data(
            ct_mean,
            x_mean,
            y_mean,
            xerr_mean,
            yerr_mean,
            outpath+'.csv',
            old_dat_path = old_dat_path,
            new_sim_path = new_sim_path,
            old_sim_path = old_sim_path,
            count_key = count_key,
            yerr_key = yerr_key,
            xs_ratio = xs_ratio,
            lumi_ratio = lumi_ratio,
            tpol_factor = tpol_factor,
            tdil_factor = tdil_factor,
            yvalue = graph_yvalue,
            xvar_keys = [xvar],
            sgasym = sgasyms[sgasym_idx],
            aliases = aliases,
        )

        # Reset graph data
        ct_mean = rescaled_graph['ct_mean']
        x_mean = rescaled_graph['x_mean']
        y_mean = rescaled_graph['y_mean']
        xerr_mean = rescaled_graph['xerr_mean']
        yerr_mean = rescaled_graph['yerr_mean']
        y_min = rescaled_graph['y_min']
        y_max = rescaled_graph['y_max']
        y_std = rescaled_graph['y_std']
        ydiff_mean = rescaled_graph['ydiff_mean']
        ydiff_std = rescaled_graph['ydiff_std']
        ydiff_min = rescaled_graph['ydiff_min']
        ydiff_max = rescaled_graph['ydiff_max']
        scaling = rescaled_graph['scaling']
        acceptanceratio = rescaled_graph['acceptanceratio']

    # Set up plot
    ax1.set_xlim(*xlims)
    ax1.set_ylim(*ylims)
    ax1.set_title(title,usetex=True)
    ax1.set_xlabel(xlabel,usetex=True)
    ax1.set_ylabel(ylabel,usetex=True)

    # Plot projection variable distribution histograms
    if len(hist_paths)>0:
        plot_hists(
            ax1,
            hist_paths = hist_paths,
            hist_keys = hist_keys,
            clone_axis = hist_clone_axis,
            ylabel = hist_ylabel,
            ylims = hist_ylims,
            histtype = histtype,
            hist_colors = hist_colors,
            alpha=hist_alpha,
            linewidth=hist_linewidth,
            density=hist_density,
            log=hist_log,
            hist_labels = hist_labels,
            binlims = binlims,
            vlinestyle = vlinestyle,
            vline_hist_idx = vline_hist_idx,
            legend_loc = hist_legend_loc,
            hist_dim = hist_dim,
        )

    # Plot systematic errors
    if yerr_syst is not None:
        g1 = ax1.errorbar(x_mean,y_mean,xerr=None,yerr=yerr_syst,
                    ecolor=fill_color, elinewidth=elinewidth*20, capsize=0,
                    color=fill_color, marker='o', linestyle=linestyle, alpha=0.5,
                    linewidth=0, markersize=0,label='Systematic error')

    # Plot standard deviation of aggregated injected values
    if y_std is not None:
        fb = ax1.fill_between(x_mean, np.add(y_mean,y_std), np.add(y_mean,-y_std), alpha=0.2, label='$\\pm1\\sigma$ Band', color=fill_color)

    # Plot results
    if yerr_mean is not None:
        g2 = ax1.errorbar(x_mean,y_mean,xerr=xerr_mean if plot_xerrors else None,yerr=yerr_mean,
                            ecolor=ecolor, elinewidth=elinewidth, capsize=capsize,
                            color=sg_colors[sgasym_idx], marker='o', linestyle=linestyle,
                            linewidth=linewidth, markersize=markersize,label=sgasym_labels[sgasym_idx])

    # Add zero line
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)

    # Draw injected asymmetries
    if show_injected_asymmetries:

        # Plot injected signal asymmetries
        plot_injected_asyms(
            ax1,
            sgasyms,
            sgasym_labels,
            sg_colors,
            sgasym_idx = sgasym_idx,
            ylims = ylims,
            label_base='Injected Signal ',
            linestyle='--',
            linewidth=axlinewidth,
        )

        # Plot injected background asymmetries
        plot_injected_asyms(
            ax1,
            bgasyms,
            bgasym_labels,
            bg_colors,
            sgasym_idx = -1,#NOTE: Don't emphasize any background asymmetries for now.
            ylims = ylims,
            label_base='Injected Background ',
            linestyle='--',
            linewidth=axlinewidth,
        )

    # Add water mark
    if watermark is not None and watermark!='': plot_watermark(ax1,watermark=watermark)

    # Plot legend
    if legend_loc is not None and legend_loc!='': ax1.legend(loc=legend_loc)

    # Check whether you have graph data to save to CSV
    if ct_mean is None: return

    # Save plot data to csv
    delimiter = ","
    cols      = ["bin","count","x","y","xerr","yerr","xerrsyst","yerrsyst"] if not rescale else ["bin","count","x","y","xerr","yerr","acceptanceratio","scaling"]
    xerr_syst_shape = np.shape(xerr_syst)
    yerr_syst_shape = np.shape(yerr_syst)
    if (len(xerr_syst_shape)==2):
        cols = ["bin","count","x","y","xerr","yerr","xerrsystlow","xerrsysthigh","yerrsyst"]
    if (len(yerr_syst_shape)==2):
        cols = ["bin","count","x","y","xerr","yerr","xerrsyst","xerrsyst","yerrsystlow","yerrsysthigh"]
    if (len(xerr_syst_shape)==2 and len(yerr_syst_shape)==2):
        cols = ["bin","count","x","y","xerr","yerr","xerrsystlow","xerrsysthigh","yerrsystlow","yerrsysthigh"]
    header    = delimiter.join(cols) #NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    fmt       = ["%.3g" for i in range(len(cols)-2)]
    fmt       = ["%d","%d",*fmt]
    comments  = ""

    # Save plot data
    save_graph_to_csv(
            outpath+'.csv' if not rescale else outpath+'_rescaled.csv',
            ct_mean,
            x_mean,
            y_mean,
            xerr=xerr_mean,
            yerr=yerr_mean,
            xerr_syst=xerr_syst if not rescale else acceptanceratio,
            yerr_syst=yerr_syst if not rescale else scaling,
            delimiter=delimiter,
            header=header,
            fmt=fmt,
            comments=comments
        )

    # Save ydiffs for MC asym injection systematics
    if ydiff_mean is not None:
        save_graph_to_csv(
            outpath+'_ydiff.csv' if not rescale else outpath+'_rescaled_ydiff.csv',
            ct_mean,
            x_mean,
            ydiff_mean,
            xerr=xerr_mean,
            yerr=ydiff_std,
            delimiter=delimiter,
            header=header,
            fmt=fmt,
            comments=comments
        )

def plot_results_array(
        graph_array,
        plot_results_kwargs_array,
        plot_results_kwargs_base = {},
        figsize = (16,10),
        outpath = 'plot_projections.pdf',
        use_grid_titles = True,
        use_grid_xlabels = True,
        use_grid_ylabels = True,
        use_grid_hist_ylabels = True,
        use_default_plt_settings = True,
    ):
    """
    Parameters
    ----------
    graph_array : list, required
        List array of graph dictionaries from :meth:`get_aggregate_graph` with the desired shape
    plot_results_kwargs_array : list, required
        List array of :meth:`plot_results` key word arguments for each graph
    plot_results_kwargs_base : dict, optional
        Dictionary of base :meth:`plot_results` key word arguments to apply to every graph
    figsize : tuple, optional
        Figure size
    outpath : str, optional
        Output graphic path
    use_grid_titles : bool, optional
        Option to assume grid titles and only use titles for plots on top row of array
    use_grid_xlabels : bool, optional
        Option to assume grid x-axis labels and only use x-axis labels for plots on bottom row of array
    use_grid_ylabels : bool, optional
        Option to assume grid y-axis labels and only use y-axis labels for plots on leftmost column of (2D) array
    use_grid_hist_ylabels : bool, optional
        Option to assume grid histogram y-axis labels and only use histogram y-axis labels for plots on rightmost column of (2D) array
    use_default_plt_settings : bool, optional
        Option to use default font and tick parameter style settings

    Description
    -----------
    Plot an array of asymmetry graph results using the :meth:`plot_results` method.
    Note that :obj:`plot_projection_kwargs` should have the same shape as :obj:`graph_array`
    and that plot and axis titles will not be set unless they are on the outside
    edge of the plot grid.  Also, note that changing the figure size to allow for 
    a :math:`(16,10)` space for each figure will be ideal when using the :obj:`use_default_plt_settings`
    option.
    """

    # Use default plotting settings
    if use_default_plt_settings: set_default_plt_settings()

    # Get and check graph matrix shape
    shape = np.shape(graph_array)
    if len(shape) not in (1,2): raise TypeError('`plot_projections()` : `graph_array` shape must have shape with len(shape) in (1,2) but shape = ',shape)

    # Create figure and axes
    f, ax = plt.subplots(*shape,figsize=figsize,squeeze=not len(shape)>1) #NOTE: squeeze will squeeze out dimension one axes!

    # Loop axes and plot results for 1D and 2D cases
    if len(shape)==1:
        for i in range(shape[0]):

                # Check for masked entry
                if graph_array[i] is None:
                    continue

                # Format graph titles and axes depending on location in grid array
                if i!=0 and use_grid_titles: plot_results_kwargs_array[i]['title'] = ''
                if i!=shape[0]-1 and use_grid_xlabels: plot_results_kwargs_array[i]['xlabel'] = ''

                # Plot results
                plot_results_kwargs = dict(plot_results_kwargs_base,**plot_results_kwargs_array[i])
                outpath_i = outpath.split('.')
                ext = outpath_i.pop(-1)
                outpath_i = '.'.join(outpath_i)
                outpath_i = ''.join([outpath_i,f'___arrloc_{i}.',ext])
                plot_results_kwargs['outpath'] = outpath_i
                plot_results(ax[i],**graph_array[i],**plot_results_kwargs)
    else:
        for i in range(shape[0]):
            for j in range(shape[1]):

                # Check for masked entry
                if graph_array[i][j] is None:
                    continue

                # Format graph titles and axes depending on location in grid array
                if j!=0 and use_grid_ylabels: plot_results_kwargs_array[i][j]['ylabel'] = ''
                if i!=0 and use_grid_titles: plot_results_kwargs_array[i][j]['title'] = ''
                if i!=shape[0]-1 and use_grid_xlabels: plot_results_kwargs_array[i][j]['xlabel'] = ''
                if j!=shape[1]-1 and use_grid_hist_ylabels: plot_results_kwargs_array[i][j]['hist_ylabel'] = ''

                # Plot results
                plot_results_kwargs = dict(plot_results_kwargs_base,**plot_results_kwargs_array[i][j])
                outpath_i = outpath.split('.')
                ext = outpath_i.pop(-1)
                outpath_i = '.'.join(outpath_i)
                outpath_i = ''.join([outpath_i,f'___arrloc_{i}_{j}.',ext])
                plot_results_kwargs['outpath'] = outpath_i
                plot_results(ax[i,j],**graph_array[i][j],**plot_results_kwargs)

    # Save figure
    rescale = plot_results_kwargs_base['lumi_ratio']>0.0 if 'lumi_ratio' in plot_results_kwargs_base else False
    if rescale: outpath = outpath.replace('.pdf','_rescaled.pdf')
    f.savefig(outpath)
