#----------------------------------------------------------------------#
# Python module for aggregating output from slurm jobs for all possible
# option combinations supplied in yaml files.
# Authors: M. McEneaney (2024, Duke University)
#----------------------------------------------------------------------#
import uproot as ur
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn

import subprocess
import os
import shutil
import yaml
import sys

def get_config_str(
        config,
        sep='_'
    ):
    """
    Parameters
    ----------
    config : required, dictionary
        Configuration mapping option names to option values
    sep : optional, string
        String separator
        Default : '_'

    Returns
    -------
    string
        String representation of configuration

    Description
    -----------
    Create a string representation of a configuration.
    """

    return (sep+sep).join([sep.join([key,sep.join([str(ele) for ele in config[key]]) if type(config[key])==list else str(config[key]) ]) for key in sorted(config)])

def get_config_out_path(
        base_dir,
        aggregate_keys,
        result_name,
        config,
        sep='_',
        ext='.pdf',
    ):
    """
    Parameters
    ----------
    base_dir : string, required
        Directory name to prepend to output file names
    aggregate_keys : list, required
        List of keys over which configurations are grouped
    result_name : string, required
        Unique string identifier for result
    config : dictionary, required
        Job configuration map
    sep : optional, string
        String separator
        Default : '_'
    sep : optional, string
        File extension
        Default : '.pdf'

    Returns
    -------
    String
        Unique output path for a file produced from the given configuration

    Description
    -----------
    Get a unique output file name for a given configuration.
    """

    job_config_name  = 'aggregate'+sep+sep+sep+(sep+sep).join([str(key) for key in sorted(aggregate_keys)])+(sep+sep+sep)
    job_config_name += get_config_str(config,sep=sep)
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
    base_dir : string, optional
        Base directory for file path
        Default : None
    base_name : string, optional
        Base name used to construct file name
        Default : ''
    binscheme_name : string, optional
        Name of binning scheme
        Default : 'x'
    ext : string, optional
        File extension
        Default : '.csv'

    Returns
    -------
    String
        Output file name

    Description
    -----------
    Get output file name for a generic kinematic binning scheme passed to an executable constructed as `<baseoutpath><binscheme_name><ext>`.
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
    configs : dictionary, required
        Map of configuration option names to option values
    aggregate_keys : list, optional
        List of keys over which to group configurations
        Default : []

    Description
    -----------
    Get a list of all possible option value combinations from a map of option names to possible values.
    """

    # Create map of elements of elements of configs and combine completely into each other for one list
    data_list = []
    for i, key in enumerate(configs):
        if key in aggregate_keys: continue
        if i==0:
                for el in configs[key]:
                    data_list.append({key: el})
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
        aggregate_keys=[]
    ):
    """
    Parameters
    ----------
    configs : dictionary, required
        Map of configuration option names to option values
    base_dir : string, required
        Path to directory in which to create job directories
    aggregate_keys : list, optional
        List of keys over which to group configurations
        Default : []

    Returns
    -------
    List
        List of job output directory names grouped across the values for `aggregate_keys`

    Description
    -----------
    Get a list of output directory names grouped across the keys listed in `aggregate_keys` for jobs generated
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
            config_str = get_config_str(config_list_i)
            job_dir = os.path.abspath(os.path.join(base_dir,config_str))
            output_dirs.append(job_dir)
            
        # Loop aggregate keys and build file list for current binning
        for key in aggregate_keys:
            for value in configs[key]:
            
                # Set the aggregate key value
                config_list_i_val = dict(config_list_i) #NOTE: Clone dictionary so that you only edit the local one.
                config_list_i_val[key] = value

                # Get job directory
                config_str = get_config_str(config_list_i_val)
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
    path : required, string
        Path to ROOT file containing histogram
    name : optional, string
        Name of `TH1` object within the ROOT file

    Returns
    -------
    np.array
        Histogram data as a numpy array or empty list if file is not found

    Description
    -----------
    Read `TH1` histogram data from a ROOT file.  Note that this will work for any histogram: (`TH1`, `TH2`, `TH3`).
    """

    # Get TH1 from ROOT file
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
    path : required, string
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
        new_path=None
    ):
    """
    Parameters
    ----------
    path : required, string
        Path to CSV file containing table
    old_path : optional, string
        Directory name to replace
        Default : None'
    new_path : optional, string
        Directory name to insert if not None
        Default : None

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
    return pd.read_csv(inpath)

def setNestedBinCuts(
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
        Default : []
    old_cut_titles : list, optional
        List of cut titles from previous recursion
        Default : []
    old_ids : list, optional
        List of bin ids from previous recursion
        Default : []
    var : string, optional
        Bin scheme variable at current recursion depth
        Default : ""
    lims_key : optional, string
        Key for bin limits
        Default : "lims"
    nested_key : optional, string
        Key for nested bin scheme
        Default : "nested"
    binvar_titles : optional, dictionary
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
        varcut_titles = [f"${varlims[idx][0]} \leq {var_title} < {varlims[idx][1]}$" for idx in range(nbins)]

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
                    setNestedBinCuts(
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

def getSchemeVars(
        node,
        nested_key = "nested"
    ):
    """
    Parameters
    ----------
    node : required, dict
        Map of bin scheme variables to bin limits arrays with either a nested or grid structure
    nested_key : optional, string
        Key for nested bin scheme structure
        Default : "nested"

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

def getNestedSchemeShape(
        node,
        lims_key   = "lims",
        nested_key = "nested"
    ):
    """
    Parameters
    ----------
    node : required, dict
        Map of bin scheme variables to bin limits arrays with a nested structure
    lims_key : optional, string
        Key for bin limits arrays
        Default : "lims"
    nested_key : optional, string
        Key for nested bin scheme structure
        Default : "nested"

    Raises
    ------
    ValueErrror:
        Raise an error if the `node` does not define a 2D nested bin structure.

    Returns
    -------
    binscheme_shape
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
    binscheme : required, dict
        Dictionary mapping binning variables to bin limits arrays
    start_idx : optional, int
        Starting integer for enumerating bin unique integer ids
        Default : 0
    id_key : optional, string
        Column name for bin unique integer ids
        Default : 'bin_id'
    binvar_titles : optional, list
        List of latex format titles for bin variables for forming cut titles
        Default : []

    Returns
    -------
    (dict, dict, pandas.DataFrame, list or None)
        Dictionary of bin ids to cuts, dictionary of bin ids to bin cut
        title dictionaries, a pandas dataframe with bin ids under
        `id_key` and the projection bin ids under the respective variable names,
        and the shape of the nested bin scheme at depth 2.  Note that the nested
        grid shape will be `None` if no 2D nested bin scheme is defined.

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
    binscheme_vars = getSchemeVars(binscheme,nested_key='nested')
    if binvar_titles is None or len(binvar_titles)!=len(binscheme_vars): binvar_titles = [var for var in binscheme_vars]

    # Check for nested bin scheme
    if (type(binscheme)==dict and 'nested' in binscheme.keys() and type(binscheme['nested'])==list and type(binscheme['nested'][0])==dict):

        # Recursively set bin scheme cuts map
        setNestedBinCuts(
            cuts,
            cut_titles,
            ids,
            binscheme,
            lims_key      = 'lims',
            nested_key    = 'nested',
            binvar_titles = binvar_titles
        )

        nested_grid_shape = getNestedSchemeShape(
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
            varcut_titles = [f"${varlims[idx][0]} \leq {var_title} < {varlims[idx][1]}$" for idx in range(nbins)]

            # Expand bin cuts and cut titles and projection ids maps
            cuts = [f"{varcut} && {cut}" for varcut in varcuts for cut in cuts] if len(cuts)>0 else varcuts
            cut_titles = [f"{var} : {varcut_title} && {cut_title}" for varcut_title in varcut_titles for cut_title in cut_titles] if len(cut_titles)>0 else [f"{var} : {varcut_title}" for varcut_title in varcut_titles]
            ids  = [dict({var:varid},**id_map) for varid in varids for id_map in ids] if len(ids)>0 else [{var:el} for el in varids]

    # Turn cuts and cuts_titles into maps
    cuts = {start_idx+idx: cuts[idx] for idx in range(len(cuts))}
    cut_titles = {start_idx+idx: {el.split(" : ")[0]:el.split(" : ")[1] for el in cut_titles[idx].split(" && ")} for idx in range(len(cut_titles))}

    # Set up data frame
    df = {id_key:[]}
    for var in binscheme_keys:
        df[var] = []

    # Add data frame entries
    for idx in range(len(ids)):
        binscheme_binid = idx + start_idx
        df[id_key].append(binscheme_binid)
        for var in binscheme_keys:
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
    grid : required, list
        List of projection ids to reshape
    nested_grid_shape : required, list
        List of (irregular) nested grid array dimensions
    fill_value : optional, int
        Fill value for padding
        Default : None

    Returns
    -------
    A reshaped grid array padded to dimension `(len(nested_grid_shape),max(nested_grid_shape))`

    Raises
    ------
    ValueError
        If `nested_grid_shape` is empty or not a list of integers or if `sum(nested_grid_shape)`
        is not the same length as the number of projections in `proj_ids`

    Description
    -----------
    Reshape projection ids for an irregular nested bin scheme padding to the largest
    nested dimension with `fill_value`.
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
        arr_var_bins={}
    ):
    """
    Parameters
    ----------
    df : required, pandas.DataFrame
        Pandas dataframe with bin ids under `id_key` and the projection bin ids under the respective variable names
    proj_vars : required, list
        Projection variables
    arr_vars : optional, list
        Variables in which to construct a grid of results
        Default : []
    id_key : optional, string
        Column name for bin unique integer ids
        Default : 'bin_id'
    arr_var_bins : optional, dict
        Dictionary of array binning scheme variable names to the specific projection bin ids to look for in those variables
        Default : {}
    nested_grid_shape : optional, list
        One dimensional list of nested grid dimensions, note this will be padded with `None` to the largest nested dimension
        Default : None

    Returns
    -------
    (all_proj_ids, arr_vars, all_proj_arr_var_ids)
        An array of each projection's unique bin ids, a list of the array binning variables encountered,
        and an array of the array variables' bin ids for each projection.  Note the array shapes should be `(*N_BINS_PROJ_VARS,*N_BINS_ARR_VARS)`.
        If `N_BINS_PROJ_VARS=[5]` and `N_BINS_ARR_VARS=[3,8]` then the shape of `all_proj_ids` and `all_proj_arr_var_ids`
        will  be `(3,8,5)`.

    Raises
    ------
    TypeError
        Raise an error if array binning variables found in keys of `arr_var_bins` are also
        found in `proj_vars` since should not simultaneously project and select a single bin.

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
    df : required, pandas.DataFrame
        Pandas dataframe containing bin ids and data
    bin_ids : required, list
        List of unique integer bin ids
    id_key : optional, string
        String identifier for bin id column
        Default : 'bin_id'
    count_key : optional, string
        String identifier for bin count column
        Default : 'count'
    xvar_keys : list, optional
        List of binning variables for which to return mean values
        Default : ['x']
    asym_key : string, optional
        Asymmetry variables for which to return mean value
        Default : 'a0'
    err_ext : string, optional
        Extension for forming error column names
        Default : 'err'

    Returns
    -------
    np.array
        Numpy array containing graph data with dimensions `(1+2*(1+N_XVAR_KEYS),*SHAPE(BIN_IDS))`

    Description
    -----------
    Read graph data (`count`,`y`,`yerr`,`x0`,`x0err`,`x1`,`x1_err`,...) for a projection plot from a pandas dataframe.
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
    graph_list : required, list
        List of graphs with dimension `(N_GRAPHS, 2*(1+N_XVAR_KEYS), N_BIN_IDS)`
    xvar_keys : list, optional
        List of binning variables for which to return mean values
        Default : ['x']
    sgasym : optional, float
        Injected signal asymmetry for computing difference of measured and injected values
        Default : 0.0

    Returns
    -------
    Dictionary
        Dictionary of mean asymmetry means and errors and other statistics names as well as bin variable means and errors to an array of their values in each kinematic bin

    Description
    -----------
    Compute the mean bin variables and asymmetry means and errors and other statistical information across a list of graphs' data from `get_graph_data()`.
    Note that in the case of that the graph dimension is greater than 1, the bin variable statistics will be returned as a list of arrays in the same order as `xvar_keys`.
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
    graph['y_std']      = np.std(graph_list[y_idx],axis=0)
    graph['y_min']      = np.min(graph_list[y_idx],axis=0)
    graph['y_max']      = np.max(graph_list[y_idx],axis=0)
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
    dfs : required, list
        Array of pandas dataframes containing bin ids and data
    proj_ids : required, list
        Array of unique integer bin ids, note if entries are set to None the graph array entry will also be set to None to allow for masked grids
    id_key : optional, string
        String identifier for bin id column
        Default : 'bin_id'
    count_key : optional, string
        String identifier for bin count column
        Default : 'count'
    xvar_keys : list, optional
        List of binning variables for which to return mean values
        Default : ['x']
    asym_key : string, optional
        Asymmetry variables for which to return mean value
        Default : 'a0'
    err_ext : string, optional
        Extension for forming error column names
        Default : 'err'
    sgasym : float or list, optional
        Injected signal asymmetry for computing difference of measured and injected values
        Default : 0.0

    Returns
    -------
    list
        An array of aggregated graph data dictionaries obtained with the same grid structure as `proj_ids`.

    Raises
    ------
    TypeError
        Raise an error the shape of the array of bin indices is not in (1,2).

    Description
    -----------
    Aggregate a set of graphs given the array of dataframes and a grid array of projection bin indices.
    Note that the returned array will have the same shape as the given `proj_ids`.
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
        xvar_keys = ['x'],
        sgasym = 0.0,
    ):
    """
    Parameters
    ----------
    ct_mean : list, required
        Count mean values for each bin
        Default : None
    x_mean : list, required
        x mean values for each bin
        Default : None
    y_mean : list, required
        y mean values for each bin
        Default : None
    xerr_mean : list, required
        x error alues for each bin
        Default : None
    yerr_mean : list, required
        y error values for each bin
        Default : None
    path : required, string
        Path where the given graph data will be stored
    old_dat_path : optional string
        Part of `path` to replace
        Default : 'old_dat_path.csv'
    new_sim_path : optional, string
        Path with which to replace `old_dat_path` to find new simulation graph CSVs
        Default : 'new_sim_path.csv'
    old_sim_path : optional, string
        Path with which to replace `old_dat_path` to find new simulation graph CSVs
        Default : 'old_sim_path.csv'
    count_key : optional, string
        CSV column key for graph bin counts
        Default : 'count'
    yerr_key : optional, string
        CSV column key for graph bin y errors
        Default : ''
    xs_ratio : optional, string
        Cross-section ratio (new/old) for scaling
        Default : 0.0
    lumi_ratio : optional, string
        Luminosity ratio (new/old) for scaling
        Default : 0.0
    xvar_keys : list, optional
        List of binning variables for which to return mean values
        Default : ['x']
    sgasym : float or list, optional
        Injected signal asymmetry for computing difference of measured and injected values
        Default : 0.0

    Returns
    -------
    dict
        A dictionary of graph data produced by `get_aggregate_graph()`

    Description
    -----------
    Rescale a graph from data (`old_dat`) loading new and old simulation graphs (`new_sim` and `old_sim`)
    from file to compute the bin dependent acceptance ratio.  Start from from either the ratio (`new_sim`/`old_sim`)
    of counts in that bin or the ratio counts estimated from asymmetry errors loaded from `yerr_key`and assuming poissonian statistics.
    The rescaling ratio is then computed by multiplying by `lumi_ratio / xs_ratio`.
    """

    # Load other graphs from csv
    new_sim_graph = load_csv(path,old_path=old_dat_path,new_path=new_sim_path)
    old_sim_graph = load_csv(path,old_path=old_dat_path,new_path=old_sim_path)

    # Get counts OR y errors from csv
    new_sim_graph_count = new_sim_graph[count_key] if yerr_key is None or yerr_key == '' else 1.0/np.square(new_sim_graph[yerr_key])
    old_sim_graph_count = old_sim_graph[count_key] if yerr_key is None or yerr_key == '' else 1.0/np.square(old_sim_graph[yerr_key])

    # Compute scaled quantities
    acceptanceratio  = np.divide(new_sim_graph_count,old_sim_graph_count)
    scaling          = acceptanceratio * lumi_ratio / xs_ratio
    scaled_ct_mean   = np.multiply(scaling,ct_mean)
    err_scaling      = np.sqrt(np.divide(ct_mean,scaled_ct_mean)) #NOTE: SCALE ERRORS ASSUMING POISSONIAN STATISTICS -> d ~ 1/sqrt(N) -> MULTIPLY BY sqrt(N_old_data/N_new_data)
    scaled_yerr_mean = np.multiply(err_scaling,yerr_mean)

    # Create a length 1 list of graph data with scaled graph results
    graph_list = np.array([[scaled_ct_mean,y_mean,scaled_yerr_mean,x_mean,xerr_mean]])

    graph = get_aggregate_graph(
            graph_list,
            xvar_keys=xvar_keys,
            sgasym=sgasym
        )

    # Set systematic errors to scaling fractions
    graph['scaling'] = scaling
    graph['acceptanceratio'] = acceptanceratio

    return graph

def get_cut_array(
        cut_titles,
        proj_ids,
        arr_vars,
    ):
    """
    Parameters
    ----------
    cuts : required, dict
        Dictionary of bin ids to list of bin cuts by variable
    proj_ids : required, list
        Array of unique integer bin ids
    arr_vars : required, list
        Variables in which to construct a grid of cuts

    Returns
    -------
    list
        An array of bin cut dictionaries with the same grid structure as `proj_ids`.

    Raises
    ------
    TypeError
        Raise an error the shape of the array of bin indices is not in (1,2).

    Description
    -----------
    Create a grid of dictionaries of array variables to bin cut titles given a
    dictionary of bin ids to cut title lists and a grid array of projection bin indices.
    Note that the returned array will have the same shape as the given `proj_ids`.
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
    args_array : required, list
        Array of argument dictionaries
    cut_array : required, List
        List array of dictionaries of array variables to bin cut titles
    proj_ids : required, list
        Array of unique integer bin ids
    arr_vars : required, list
        Variables in which to construct a grid of cuts

    Returns
    -------
    list
        An array of bin cut dictionaries with the same grid structure as `proj_ids`.

    Raises
    ------
    TypeError
        Raise an error the shape of the arrays do not match or the length of the shapes is not in (1,2).

    Description
    -----------
    Add bin cut 'title' and 'ylabel' arguments to each argument dictionary in a grid array
    from a cut array produced from `get_cut_array`.
    Note that the returned array will have the same shape as the given `args_array`.
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
    df : required, pandas.DataFrame
        Pandas dataframe of bin migration matrix
    id_gen_key : optional, string
        Key for generated bin indices
        Default : 'binid_gen'
    id_rec_key : optional, string
        Key for reconstructed bin indices
        Default : 'binid_rec'
    mig_key : optional, string
        Key for bin migration fractions
        Default : 'mig'

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
    df : required, pandas.DataFrame
        Pandas dataframe containing bin ids and data

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe subset containing only elements whose `id_key` entries are in `bin_ids`

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
        Array describing a 1D graph with structure `graph[{x_mean,x_err,y_mean,y_err,...},{nbins}]`
    offset : float, required
        Value by which to offset graph values
    axis : int, optional
        Axis along which to offset the graph
        Default : 0

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
    filename : required, string
        Output file name
    data : array, required
        2D data array with dimensions `[N_COLUMNS,N_ROWS]` TODO CHECK THIS!!!
    delimiter : optional, string
        CSV format delimiter
        Default : ","
    header : optional, string
        CSV header
        Default : None
    fmt : optional, string
        CSV column formats
        Default : None
    comments : optional, string
        CSV comments
        Default : ""

    Description
    -----------
    Save a square data array of dimensions `[N_COLUMNS,N_ROWS]` to a text file.  TODO CHECK THIS!!!
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
    filename : required, string
        Output file name
    ct : required, list
        Graph count values
    x : required, list
        Graph x values
    y : required, list
        Graph y values
    xerr : optional, list
        Graph x error values
        Default : None
    yerr : optional, list
        Graph y error values
        Default : None
    xerr_syst : optional, list
        Graph x systematic error values
        Default : None
    yerr_syst : optional, list
        Graph y systematic error values
        Default : None
    delimiter : optional, string
        CSV format delimiter
        Default : ","
    header : optional, string
        CSV header
        Default : None
    fmt : optional, string
        CSV column formats
        Default : None
    comments : optional, string
        CSV comments
        Default : ""

    Description
    -----------
    Write a graph to a CSV file with optional errors and systematic errors.
    """

    # Create data array
    data = []
    if ct is None or len(ct)==0: ct = [0.0 for el in x]
    if xerr is None or len(xerr)==0: xerr = [0.0 for el in x]
    if yerr is None or len(yerr)==0: yerr = [0.0 for el in x]
    if xerr_syst is None or len(xerr_syst)==0: xerr_syst = [0.0 for el in x]
    if yerr_syst is None or len(yerr_syst)==0: yerr_syst = [0.0 for el in x]
    for i, el in enumerate(x):
        data.append([i, ct[i], x[i], y[i], xerr[i], yerr[i], xerr_syst[i], yerr_syst[i]])
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
    filename : required, string
        Output file name
    x : required, list
        Graph x values
    yerr_syst : optional, list
        Graph y systematic error values decomposed into the different sources of systematic error
        Default : None
    delimiter : optional, string
        CSV format delimiter
        Default : ","
    header : optional, string
        CSV header
        Default : None
    fmt : optional, string
        CSV column formats
        Default : None
    comments : optional, string
        CSV comments
        Default : ""

    Description
    -----------
    Write a set of graph y systematic errors to a CSV file with the systematic error values broken down by source.
    This means the argument `yerr_syst` should have shape (`nbins`,`nsources`) where `nbins` is the number of 
    kinematic bins and `nsources` is the number of sources of systematic error.
    """

    # Create data array
    data = []
    if yerrs_syst is None or len(yerrs_syst)==0: yerrs_syst = [[0.0] for el in x]
    for i, el in enumerate(x):
        data.append([i, x[i], *yerrs_syst[i]])
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
    bin_mig_mat : required, np.array
        2D bin migration matrix with (i,j) `->` (generated,reconstructed)
    base_dir : string, required
        Path to directory in which matrix will be saved
        Default : 'systematics/bin_migration/'
    basename : optional, string
        Name of reconstructed bin variable
        Default : 'x'
    delimiter : optional, string
        CSV format delimiter
        Default : ","
    header : optional, string
        CSV header
        Default : None
    fmt : optional, string
        CSV column formats, automatically set if not specified
        Default : None
    comments : optional, string
        CSV comments
        Default : ""

    Raises
    ------
    TypeError
        Raise an error if the bin migration matrix is not square

    Description
    -----------
    Save a 2D bin migration matrix mapping generated bins to reconstructed bins to a CSV file
    with an added initial row and column for the bin indices.  Note that files will be saved
    to `<base_dir>/bin_mig_mat_<basename>.csv`.
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
    df : required, pandas.DataFrame
        Pandas dataframe of asymmetry results
    inv_bin_mig_mat : required, np.array
        Inverse of bin migration matrix mapping generated bins to reconstructed bins
    results_keys : optional, list
        List of keys to results entries to which to apply bin migration corrections
        Default : 'a0'

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
        Default : None
    systematic_scales_mat : np.array, optional
        Array of systematic errors to be scaled by y values before being added in quadrature to other systematic errors
        Default : None
    systematics_additive_mat : np.array, optional
        Array of absolute systematic errors to add to other systematic errors
        Default : None

    Returns
    -------
    np.array
        Array of systematics values added in quadrature

    Description
    -----------
    Compute the systematic errors for a 1D binning scheme given any combination
    of bin migration matrix (this will be inverted internally with `np.linalg.inv`),
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
        Default : 0
    ylims : tuple, optional
        y limits for plotting
        Default : (-1.0,1.0)
    label_base : string, optional
        Base label to prepend to ytitles for injected asymmetries
        Default : 'Injected Signal '
    linestyle : string, optional
        Line style
        Default : '--'
    linewidth : int, optional
        Line width
        Default : 1

    Description
    -----------
    Plot the injected asymmetries for each bin in a 1D binning scheme offsetting repeat values
    by a small amount.  Note that injected asymmetries should have shape `(N_ASYMS)` if plotting
    constant asymmetries or `(N_ASYMS,2,N_POINTS)` if you would like to plot function data `(x,y)`.
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
    watermark : string, optional
        Watermark text
        Default : 'CLAS12 Preliminary'

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
    linestyle : string, optional
        Line style
        Default : 'dotted'

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
        hist_labels = None,
        binlims = [],
        vlinestyle = 'dotted',
        vline_hist_idx = -1,
        legend_loc = 'upper right',
    ):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    clone_axis : Boolean, optional
        Option to create a twin y axis sharing the x axis of the given axis
        Default : True
    ylabel : string, optional
        Y axis label for twin axis
        Default : 'Density'
    ylims : tuple, optional
        Y axis limits for twin axis
        Default : (0.0,0.05)
    histtype : string, optional
        Matplotlib.pyplot histogram type
        Default : 'step'
    hist_colors : list, optional
        List of histogram colors
        Default : None
    alpha : float, optional
        Alpha plotting paramater for histograms
        Default : 0.5
    linewidth : int, optional
        Line width for plotting histograms
        Default : 2
    density : boolean, optional
        Option to normalize histograms
        Default : True
    hist_labels : list, optional
        List of histogram labels
        Default : None
    binlims : list, optional
        List of bin limits in a 1D binning scheme
        Default : []
    vlinestyle : string, optional
        Vertical line style for drawing bin limits on histogram
        Default : 'dotted'
    vline_hist_idx : int, optional
        Index of histogram for which to draw vertical lines for bin limits
        Default : -1
    legend_loc : string, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to None or ''
        Default : 'upper right'

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
            density=False
        )

        # Plot bin limits if supplied and we are on the last histogram
        if idx==(vline_hist_idx if vline_hist_idx>=0 else len(hist_paths)-1) and len(binlims)>0:
            plot_vlines(
                h,
                binlims,
                linestyle = vlinestyle,
            )

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
        Default : ['x','Q2']
    cols : list, optional
        List of kinematics LaTeX titles
        Default : ['x','Q^{2}']
    err_ext : string, optional
        Extension for forming column names of kinematic variable errors
        Default : '_err'
    sep : string, optional
        Separator string for kinematics values in bin title
        Default : ' , '

    Returns
    -------
    string
        Title string

    Description
    -----------
    Create a title string for the requested bin index and kinematic variables showing the mean and error values for
    each of the requested kinematic variables in the given bin.
    """

    return sep.join(
        [
            f"$<{col_titles[col]}> = {df.iloc[bin_id].loc[col]:.2f}\pm{df.iloc[bin_id].loc[cols[idx]+err_ext]:.2f}$"
            for idx, col in enumerate(cols)
        ]
    )

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
        ytitle  = '$\Delta \mathcal{A}$',
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
        figsize = (16,10),
    ):
    """
    Parameters
    ----------
    x_means : list, required
        Mean x values for each bin
    yerr_syst : np.array, required
        Array of absolute systematic error in each bin further indexed by the sources of systematic error
    palette : string, optional
        Seaborn color palette
        Default : 'Dark2'
    stacked : Boolean, optional
        Whether to stack histograms from different sources of systematic error
        Default : False
    syst_names : List, optional
        List of column names for each source of systematic error
        Default : None
    syst_labels : List, optional
        List of labels for each source of systematic error
        Default : None
    xlims : tuple, optional
        x limits for plotting
        Default : (0.0,1.0)
    ylims : tuple, optional
        y limits for plotting
        Default : (-1.0,1.0)
    xvar : string, optional
        Bin variable name
        Default : 'x'
    title : string, optional
        Plot title
        Default : 'Systematic Errors'
    xtitle : string, optional
        x axis title
        Default : '$Q^{2} (GeV^{2})$'
    ytitle : string, optional
        y axis title
        Default : '$\Delta \mathcal{A}$'
    outpath : string, optional
        Name of output pdf
        Default : 'systematics.pdf'
    watermark : String, optional
        Optional watermark to put on top of plot
        Default : 'CLAS12 Preliminary'
    use_default_plt_settings : Boolean, optional
        Option to use default font and tick parameter style settings
        Default : True
    legend_loc : string, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to None or ''
        Default : 'upper left'
    ecolor : string, optional
        Error line color
        Default : 'black'
    ecolor : float, optional
        Error line width
        Default : 2.0
    capsize : int, optional
        Error cap size
        Default : 18
    capthick : float, optional
        Error cap thickness
        Default : 2.0
    marker : string, optional
        Marker type
        Default : 'o'
    markersize : int, optional
        Marker size
        Default : 20
    linestyle : string, optional
        Line style
        Default : None
    linewidth : float, optional
        Line width
        Default : 0.0
    gridlinewidth : float, optional
        Grid line width
        Default : 0.5
    axlinewidth : float, optional
        Axis line and injected asymmetries line width
        Default : 1.0
    figsize : tuple, optional
        Figure size
        Default : (16,10)

    Description
    -----------
    Plot the systematic error for each bin in a 1D binning scheme broken down by sources of systematic error.
    Save systematics breakdowns to CSV in `<outpath>.csv`.
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
    s1 = plt.hist(xbins, weights=yerr_syst, bins=nbins, alpha=0.5, label=syst_labels, stacked=stacked)

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
        ylabel = '$\mathcal{A}$',
        sgasym_labels = ['$\mathcal{A}$'],
        bgasym_labels = ['$\mathcal{A}$'],
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
        hist_labels = None,
        binlims = [],
        vlinestyle = 'dotted',
        vline_hist_idx = -1,
        hist_legend_loc = 'upper right',
        old_dat_path = 'old_dat_path.csv',
        new_sim_path = 'new_sim_path.csv',
        old_sim_path = 'old_sim_path.csv',
        count_key = 'count',
        yerr_key = '',
        xs_ratio = 1.0,
        lumi_ratio = 0.0,
    ):
    """
    Parameters
    ----------
    ax1 : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot figure axis
    ct_mean : list, optional
        Count mean values for each bin
        Default : None
    x_mean : list, optional
        x mean values for each bin
        Default : None
    y_mean : list, optional
        y mean values for each bin
        Default : None
    xerr_mean : list, optional
        x error alues for each bin
        Default : None
    yerr_mean : list, optional
        y error values for each bin
        Default : None
    xerr_syst : list, optional
        x systematic error alues for each bin
        Default : None
    yerr_syst : list, optional
        y systematic error values for each bin
        Default : None
    y_min : list, optional
        y minimum values for each bin
        Default : None
    y_max : list, optional
        y maximum values for each bin
        Default : None
    y_std : list, optional
        y standard deviation values for each bin
        Default : None
    ydiff_mean : list, optional
        y difference from injected signal asymmetry mean values for each bin
        Default : None
    ydiff_std : list, optional
        y difference from injected signal asymmetry standard deviation values for each bin
        Default : None
    ydiff_min : list, optional
        y difference from injected signal asymmetry minimum values for each bin
        Default : None
    ydiff_max : list, optional
        y difference from injected signal asymmetry maximum values for each bin
        Default : None
    xlims : tuple, optional
        x limits for plotting
        Default : (0.0,1.0)
    ylims : tuple, optional
        y limits for plotting
        Default : (-1.0,1.0)
    title : string, optional
        Plot title
        Default : 'x'
    xvar : string, optional
        Bin variable name
        Default : 'x'
    xlabel : string, optional
        x axis label
        Default : '$x$'
    ylabel : string, optional
        y axis label
        Default : '$\mathcal{A}$'
    sgasym_labels : list, optional
        List of signal asymmetry labels
        Default : ['$\mathcal{A}$']
    bgasym_labels : list, optional
        List of background asymmetry labels
        Default : ['$\mathcal{A}$']
    sgasym_idx : int, optional
        Index of injected signal asymmetry
        Default : 0
    sgasyms : list, optional
        List of injected signal asymmetries
        Default : [0.10]
    bgasyms : list, optional
        List of injected background asymmetries
        Default : [0.00]
    sg_colors : list, optional
        List of signal asymmetry plotting colors
        Default : ['blue']
    bg_colors : list, optional
        List of background asymmetry plotting colors
        Default : ['red']
    fill_color : string, optional
        Color of 1 sigma band or systematic uncertainties
        Default : 'gray'
    outpath : string, optional
        Name of output pdf
        Default : 'systematics.pdf'
    watermark : String, optional
        Optional watermark to put on top of plot
        Default : 'CLAS12 Preliminary'
    show_injected_asymmetries : Boolean, optional
        Option to show injected signal and background asymmetries
        Default : False
    legend_loc : string, optional
        Matplotlib.pyplot legend location string, will not be plotted if set to None or ''
        Default : 'upper left'
    ecolor : string, optional
        Error line color
        Default : 'black'
    ecolor : float, optional
        Error line width
        Default : 2.0
    capsize : int, optional
        Error cap size
        Default : 18
    capthick : float, optional
        Error cap thickness
        Default : 2.0
    marker : string, optional
        Marker type
        Default : 'o'
    markersize : int, optional
        Marker size
        Default : 20
    linestyle : string, optional
        Line style
        Default : None
    linewidth : float, optional
        Line width
        Default : 0.0
    gridlinewidth : float, optional
        Grid line width
        Default : 0.5
    axlinewidth : float, optional
        Axis line and injected asymmetries line width
        Default : 1.0
    hist_clone_axis : Boolean, optional
        Option to create a twin y axis sharing the x axis of the given axis
        Default : True
    hist_ylabel : string, optional
        Y axis label for twin axis
        Default : 'Density'
    hist_ylims : tuple, optional
        Y axis limits for twin axis
        Default : (0.0,0.05)
    histtype : string, optional
        Matplotlib.pyplot histogram type
        Default : 'step'
    hist_colors : list, optional
        List of histogram colors
        Default : None
    hist_alpha : float, optional
        Alpha plotting paramater for histograms
        Default : 0.5
    hist_linewidth : int, optional
        Line width for plotting histograms
        Default : 2
    density : boolean, optional
        Option to normalize histograms
        Default : True
    hist_labels : list, optional
        List of histogram labels
        Default : None
    binlims : list, optional
        List of bin limits in a 1D binning scheme
        Default : []
    vlinestyle : string, optional
        Vertical line style for drawing bin limits on histogram
        Default : 'dotted'
    vline_hist_idx : int, optional
        Index of histogram for which to draw vertical lines for bin limits
        Default : -1
    hist_legend_loc : string, optional
        Matplotlib.pyplot legend location string for histograms, will not be plotted if set to None or ''
        Default : 'upper right'
    old_dat_path : optional string
        Part of `path` to replace
        Default : 'old_dat_path.csv'
    new_sim_path : optional, string
        Path with which to replace `old_dat_path` to find new simulation graph CSVs
        Default : 'new_sim_path.csv'
    old_sim_path : optional, string
        Path with which to replace `old_dat_path` to find new simulation graph CSVs
        Default : 'old_sim_path.csv'
    count_key : optional, string
        CSV column key for graph bin counts
        Default : 'count'
    yerr_key : optional, string
        CSV column key for graph bin y errors
        Default : ''
    xs_ratio : optional, string
        Cross-section ratio (new/old) for scaling
        Default : 0.0
    lumi_ratio : optional, string
        Luminosity ratio (new/old) for scaling
        Default : 0.0

    Description
    -----------
    Plot asymmetry results for each bin in a 1D binning scheme showing projection variable histograms, bin limits, systematic errors,
    a standard deviation band for aggregate graphs, injected asymmetries, watermark, and legend if desired.  Save results
    and differences from the injected signal to CSV in `<outpath>.csv` and `<outpath>_ydiff.csv`.
    Optionally, rescale the input graph using `rescale_graph_data`.
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
            xvar_keys = [xvar],
            sgasym = sgasyms[sgasym_idx],
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
            hist_labels = hist_labels,
            binlims = binlims,
            vlinestyle = vlinestyle,
            vline_hist_idx = vline_hist_idx,
            legend_loc = hist_legend_loc,
        )

    # Plot systematic errors
    if yerr_syst is not None:
        g1 = ax1.errorbar(x_mean,y_mean,xerr=None,yerr=yerr_syst,
                    ecolor=fill_color, elinewidth=elinewidth*20, capsize=0,
                    color=fill_color, marker='o', linestyle=linestyle, alpha=0.5,
                    linewidth=0, markersize=0,label='Systematic error')

    # Plot standard deviation of aggregated injected values
    if y_std is not None:
        fb = ax1.fill_between(x_mean, np.add(y_mean,y_std), np.add(y_mean,-y_std), alpha=0.2, label='$\pm1\sigma$ Band', color=fill_color)

    # Plot results
    if yerr_mean is not None:
        g2 = ax1.errorbar(x_mean,y_mean,xerr=xerr_mean,yerr=yerr_mean,
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
    header    = delimiter.join(cols) #NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    fmt       = ["%.3g" for i in range(len(cols)-2)]
    fmt       = ["%d","%d",*fmt]
    comments  = ""

    # Save plot data
    save_graph_to_csv(
            outpath+'.csv',
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
            outpath+'_ydiff.csv',
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
        List array of graph dictionaries from `get_aggregate_graph()` with the desired shape
    plot_results_kwargs_array : list, required
        List array of `plot_results()` key word arguments for each graph
    plot_results_kwargs_base : dict, optiional
        Dictionary of base `plot_results()` key word arguments to apply to every graph
        Default : {}
    figsize : tuple, optional
        Figure size
        Default : (16,10)
    outpath : string, optional
        Output graphic path
        Default : 'plot_projections.pdf'
    use_grid_titles : Boolean, optional
        Option to assume grid titles and only use titles for plots on top row of array
        Default : True
    use_grid_xlabels : Boolean, optional
        Option to assume grid x-axis labels and only use x-axis labels for plots on bottom row of array
        Default : True
    use_grid_ylabels : Boolean, optional
        Option to assume grid y-axis labels and only use y-axis labels for plots on leftmost column of (2D) array
        Default : True
    use_grid_hist_ylabels : Boolean, optional
        Option to assume grid histogram y-axis labels and only use histogram y-axis labels for plots on rightmost column of (2D) array
        Default : True
    use_default_plt_settings : Boolean, optional
        Option to use default font and tick parameter style settings
        Default : True

    Description
    -----------
    Plot an array of asymmetry graph results using the `plot_results()` method.
    Note that `plot_projection_kwargs` should have the same shape as `graph_array`
    and that plot and axis titles will not be set unless they are on the outside
    edge of the plot grid.  Also, note that changing the figure size to allow for 
    a (16,10) space for each figure will be ideal when using the `use_default_plt_settings`
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
    f.savefig(outpath)
