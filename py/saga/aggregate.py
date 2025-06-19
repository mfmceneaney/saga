"""
This module implements methods for aggregating output
from slurm jobs for all possible option combinations
supplied in yaml files. It also offers methods for 
manipulating, plotting, and saving the outputs.

# Author: Matthew F. McEneaney (2024, Duke University)
"""
import os
import numpy as np
import pandas as pd
from .data import load_csv

def get_config_str(
        config,
        sep='_',
        aliases=None
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

    # Check arguments
    if aliases is None:
        aliases = {}

    return (sep+sep).join([
                aliases[key][config[key]] if (key in aliases and (type(config[key]) in (str,float,int)) and config[key] in aliases[key])
                else aliases[key][str(config[key])] if (key in aliases and (type(config[key]) not in (str,float,int)) and str(config[key]) in aliases[key])
                else sep.join([
                    key,sep.join([str(ele) for ele in config[key]]) if isinstance(config[key],list) else str(config[key])
                ]) for key in sorted(config)
                ])

def get_config_out_path(
        base_dir,
        aggregate_keys,
        result_name,
        config,
        sep='_',
        aliases=None,
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
        aggregate_keys=None
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

    # Check arguments
    if aggregate_keys is None:
        aggregate_keys = {}

    # Create map of elements of elements of configs and combine completely into each other for one list
    data_list = []
    first = True
    for key in configs:
        if key in aggregate_keys:
            continue
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
        aggregate_keys=None,
        aliases=None
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

    # Check arguments
    if aggregate_keys is None:
        aggregate_keys = []

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

def set_nested_bin_cuts(
        cuts,
        cut_titles,
        ids,
        node,
        old_cuts       = None,
        old_cut_titles = None,
        old_ids        = None,
        var            = "",
        lims_key       = "lims",
        nested_key     = "nested",
        binvar_titles  = None
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

    # Check arguments
    if old_cuts is None:
        old_cuts = []
    if old_cut_titles is None:
        old_cut_titles = []
    if old_ids is None:
        old_ids = []
    if binvar_titles is None:
        binvar_titles = {}

    # Check the YAML node
    if isinstance(node,dict):

        # Check for bin limits
        lims = []
        if lims_key in node.keys() and isinstance(node[lims_key],list):
            lims = node[lims_key]

        # Set nbins lower limit to 0 since you allow passing limits with length 0
        nbins = max(len(lims)-1,0)

        # Get bin variable title
        var_title = var if var not in binvar_titles else binvar_titles[var]

        # Loop bins and get bin cuts
        varlims = [[lims[idx],lims[idx+1]] for idx in range(nbins)]
        varids  = list(range(nbins))
        varcuts = [f"({var}>={varlims[idx][0]} && {var}<{varlims[idx][1]})" for idx in range(nbins)]
        varcut_titles = [f"${varlims[idx][0]} \\leq {var_title} < {varlims[idx][1]}$" for idx in range(nbins)]

        # Expand bin cuts and cut titles and projection ids maps
        old_cuts = [f"{varcut} && {cut}" for varcut in varcuts for cut in old_cuts] if len(old_cuts)>0 else varcuts
        old_cut_titles = [f"{var} : {varcut_title} && {cut_title}" for varcut_title in varcut_titles for cut_title in old_cut_titles] \
                    if len(old_cut_titles)>0 else [f"{var} : {varcut_title}" for varcut_title in varcut_titles]
        old_ids = [dict({var:varid},**id_map) for varid in varids for id_map in old_ids] if len(old_ids)>0 else [{var:el} for el in varids]

        # Check for nested binning
        if nested_key in node and isinstance(node[nested_key],list):

            # Get nested YAML node
            node_nested = node[nested_key]

            # Loop nested bins
            for ibin, node_nested_bin in enumerate(node_nested): #NOTE: These are not bin limits just maps to bin limits for each bin, so loop normally.

                # Loop nested bin variables (only expect one!)
                for it_key in node_nested_bin:

                    # Get nested YAML node
                    it_nested = node_nested_bin[it_key]

                    # Create a new vector for uniqueness along different recursion branches
                    new_old_cuts = []
                    new_old_cut_titles = []
                    new_old_ids = []
                    if ibin<len(old_cuts):
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
            for ibin, old_cuts_i in enumerate(old_cuts):
                cuts.append(old_cuts_i)
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
    if isinstance(node,dict) and nested_key in node and isinstance(node[nested_key],list) and isinstance(node[nested_key][0],dict):

        # Get nested node
        node_nested = dict(node)

        # Loop nested nodes
        while isinstance(node_nested,dict) and nested_key in node_nested and isinstance(node_nested[nested_key],list) and isinstance(node_nested[nested_key][0],dict):

            binvar = list(node_nested[nested_key][0].keys())[0]
            binscheme_vars.append(binvar)
            node_nested = dict(node_nested[nested_key][0][binvar])

        return binscheme_vars

    # Default to case you have a grid scheme
    return list(node.keys())

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
    if not (isinstance(node,dict) and nested_key in node and isinstance(node[nested_key],list) and isinstance(node[nested_key][0],dict)):
        raise ValueError(f"`node` must have a 2D nested bin scheme structure, but node = {node}")

    # Get nested node
    node_nested = dict(node)

    # Loop nested node
    depth = 0
    while isinstance(node_nested,dict) and nested_key in node_nested and isinstance(node_nested[nested_key],list) and isinstance(node_nested[nested_key][0],dict):

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


def get_binscheme_cuts_and_ids(
        binscheme,
        start_idx=0,
        id_key='bin_id',
        binvar_titles = None,
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
    if binvar_titles is None or len(binvar_titles)!=len(binscheme_vars):
        binvar_titles = list(binscheme_vars)

    # Check for nested bin scheme
    if isinstance(binscheme,dict) and 'nested' in binscheme and isinstance(binscheme['nested'],list) and isinstance(binscheme['nested'][0],dict):

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
            if isinstance(binscheme[var],dict) and 'nbins' in binscheme[var] and 'lims' in binscheme[var]:
                nbins = binscheme[var]['nbins']
                lims = [(binscheme[var]['lims'][1]-binscheme[var]['lims'][0])/nbins * i + binscheme[var]['lims'][0] for i in range(nbins+1)]
            varlims = [[lims[idx],lims[idx+1]] for idx in range(nbins)]
            varids  = list(range(nbins))
            varcuts = [f"({var}>={varlims[idx][0]} && {var}<{varlims[idx][1]})" for idx in range(nbins)]
            varcut_titles = [f"${varlims[idx][0]} \\leq {var_title} < {varlims[idx][1]}$" for idx in range(nbins)]

            # Expand bin cuts and cut titles and projection ids maps
            cuts = [f"{varcut} && {cut}" for varcut in varcuts for cut in cuts] if len(cuts)>0 else varcuts
            cut_titles = [f"{var} : {varcut_title} && {cut_title}" for varcut_title in varcut_titles for cut_title in cut_titles] \
                        if len(cut_titles)>0 else [f"{var} : {varcut_title}" for varcut_title in varcut_titles]
            ids  = [dict({var:varid},**id_map) for varid in varids for id_map in ids] if len(ids)>0 else [{var:el} for el in varids]

    # Turn cuts and cuts_titles into maps
    cuts = {start_idx+idx: cuts[idx] for idx in range(len(cuts))}
    cut_titles = {start_idx+idx: {el.split(" : ")[0]:el.split(" : ")[1] for el in cut_titles[idx].split(" && ")} for idx in range(len(cut_titles))}

    # Set up data frame
    df = {id_key:[]}
    for var in binscheme_vars:
        df[var] = []

    # Add data frame entries
    for idx, el in enumerate(ids):
        binscheme_binid = idx + start_idx
        df[id_key].append(binscheme_binid)
        for var in binscheme_vars:
            df[var].append(el[var])

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
    if not isinstance(nested_grid_shape,list):
        raise TypeError('`nested_grid_shape` must be a list of integers')
    if len(nested_grid_shape)==0:
        raise ValueError('`nested_grid_shape` must be a non-empty list')
    if not isinstance(nested_grid_shape[0],int):
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
        arr_vars = None,
        id_key='bin_id',
        arr_var_bins=None,
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

    # Check arguments
    if arr_vars is None:
        arr_vars = []
    if arr_var_bins is None:
        arr_var_bins = {}

    # Check projection variables
    if len(proj_vars)>2:
        print('WARNING: `get_projection_ids` : Are you sure you want more than 2 projection variables?')

    # Get list of grouping variables
    for key in df.keys():
        if key!=id_key and key not in proj_vars and key not in arr_var_bins and len(arr_vars)==0:
            arr_vars.append(key)

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
    if arr_var_bins_cut!='':
        query = ' and '.join([query,arr_var_bins_cut])
    start_bins_slice = df.query(query)

    # Loop starting bins and create projection lists
    all_proj_ids = []
    all_proj_arr_var_ids = []
    for proj_idx in range(len(start_bins_slice)):

        # Find bins that match starting bin in array bin variable indices
        arr_var_ids = [start_bins_slice[arr_var].iloc[proj_idx] for arr_var in arr_vars]
        query = ' and '.join([f'{arr_vars[idx]}=={arr_var_ids[idx]}' for idx in range(len(arr_var_ids))])
        if arr_var_bins_cut!='':
            query = ' and '.join([query,arr_var_bins_cut])
        proj_ids = df.query(query)
        proj_ids = proj_ids.sort_values(proj_vars)[id_key].values #NOTE: SORT BY PROJECTION VARIABLE BINS
        proj_ids = proj_ids.reshape(nbins) #NOTE: RESIZE BY APPROPRIATE NUMBER OF BINS IF REQUESTED.
        all_proj_ids.append(proj_ids.tolist())
        all_proj_arr_var_ids.append(arr_var_ids)

    all_proj_ids = np.reshape(all_proj_ids,grid_shape) if nested_grid_shape is None else reshape_nested_grid(all_proj_ids, nested_grid_shape)
    all_proj_arr_var_ids = np.reshape(all_proj_arr_var_ids,(*nbins_arr,len(arr_vars)))

    return all_proj_ids, arr_vars, all_proj_arr_var_ids

def get_graph_data(
                    df,
                    bin_ids,
                    id_key='bin_id',
                    count_key='count',
                    xvar_keys=None,
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

    # Check arguments
    if xvar_keys is None:
        xvar_keys = []

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
        xvar_keys=None,
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

    # Check arguments
    if xvar_keys is None:
        xvar_keys = []

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
        for xvar_idx in range(len(xvar_keys)):
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
        xvar_keys=None,
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

    # Check shape
    if not len(shape) in (2,3):
        raise TypeError(f"`get_graph_array` : `proj_ids` must have len(shape) in (2,3) but shape = {shape}")

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
                sgasym=sgasym[i][j] if not isinstance(sgasym,float) else sgasym
            ) if proj_ids[i][j] is not None else None for j in range(shape[1])] for i in range(shape[0]) #NOTE Allow masked grid
        ]

    # Create a graph array in the 1D grid case
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
            sgasym=sgasym[i] if not isinstance(sgasym,float) else sgasym
        ) for i in range(shape[0])
    ]


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
        xvar_keys = None,
        sgasym = 0.0,
        aliases=None,
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
    if yvalue>=-1:
        scaled_yerr_mean *= np.sqrt(1-np.square(yvalue*tpol_factor))

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
        config = None,
        aggregate_config = None,
        chain_configs = None,
        aliases=None,
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
    new_sim_df = load_csv(path,old_path=old_dat_path,new_path=new_sim_path,aliases=aliases)
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
    if yvalue>=-1:
        scaled_yerrs *= np.sqrt(1-np.square(yvalue*tpol_factor))

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
    ):
    """
    Parameters
    ----------
    cuts : dict, required
        Dictionary of bin ids to list of bin cuts by variable
    proj_ids : list, required
        Array of unique integer bin ids

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

    # Check shape
    if len(shape) not in (2,3):
        raise TypeError('`get_cut_array` : `proj_ids` must have len(shape) in (2,3) but shape = ',shape)

    #NOTE: All BINS IN A 1D BINNING PROJECTION SHOULD HAVE THE SAME BIN CUTS FOR THE ARRAY VARIABLES.

    # Create a graph array in the 2D grid case
    if len(shape)==3:
        return [[
            cut_titles[proj_ids[i][j][0]]
            for j in range(shape[1])] for i in range(shape[0])
        ]

    # Create a graph array in the 1D grid case
    return [
        cut_titles[proj_ids[i][0]]
        for i in range(shape[0])
    ]


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
    if shape!=np.shape(cut_array):
        raise TypeError(
            "`add_cut_array` : `args_array` must have same shape as `cut_array` " +
            f"but shapes are {shape} and {np.shape(cut_array)}"
        )

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
        raise TypeError(
            "`get_bin_migration_matrix` : Generated and reconstructed unique bin ids must have same members " +
            f"but found: unique_binids_gen={unique_binids_gen}, unique_binids_rec={unique_binids_rec}"
        )

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

def apply_bin_mig(
        df,
        inv_bin_mig_mat,
        results_keys = (),
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
        bin_migration_mat_inv = np.linalg.inv(bin_migration_mat)
        new_systematics = np.add(results,-np.matmul(bin_migration_mat_inv,results)) # DeltaA = a - f_inv . a
        systematics = np.sqrt(np.square(systematics) + np.square(new_systematics))

    # Apply multiplicative scale systematics, note that these should already be summed over all sources of systematic error
    if systematic_scales_mat is not None:
        systematics = np.sqrt(np.square(systematics) + np.square(np.multiply(results,systematic_scales_mat))) #NOTE: IMPORTANT!  ADD IN QUADRATURE.

     # Apply additive scale systematics, note that these should already be summed over all sources of systematic error
    if systematics_additive_mat is not None:
        systematics += systematics_additive_mat

    return systematics
