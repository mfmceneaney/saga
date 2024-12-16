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

def get_list(divisions,aggregate_keys=[]):
    """
    Parameters
    ----------
    divisions : dictionary, required
        Map of option names to option values
    aggregate_keys : list, optional
        List of keys over which to group configurations
        Default : []

    Description
    -----------
    Get a list of all possible option value combinations from a map of option names to possible values.
    """

    # Create map of elements of elements of divisions and combine completely into each other for one list
    data_list = []
    for i, key in enumerate(divisions):
        if key in aggregate_keys: continue
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

def get_out_file_list_1D_binning(divisions,base_dir,submit_path,yaml_path,var_lims,get_out_file_name,aggregate_keys=[]):

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
    var_lims : dictionary, required
        Map of bin variable names to minimum and maximum bounds
    get_out_file_name : function, required
        Function to generate output file name from configuration map and additional 1D binning arguments `binvar`, `binvar_min`, `binvar_max`
    aggregate_keys : list, optional
        List of keys over which to group configurations
        Default : []

    Returns
    -------
    List
        List of job output file names grouped across the values for `aggregate_keys`

    Description
    -----------
    Get a list of output file names grouped across the keys listed in `aggregate_keys` for jobs generated
    from all possible option value combinations from a map of option names to possible values.
    """
    
    # Get a list of all possible option value combinations from divisions
    data_list = get_list(divisions,aggregate_keys=aggregate_keys)
    
    # Create list for output file names
    out_file_list = []

    # Loop resulting list
    for _data_list_i in data_list:

        # Add in aggregate keys
        data_list_i = dict(_data_list_i)
        
        # Loop binning variables first
        for binvar in var_lims:
            binvar_min = var_lims[binvar][0]
            binvar_max = var_lims[binvar][1]
            data_list_i_binvar = dict(_data_list_i)
            data_list_i_binvar['binvar'] = binvar
            data_list_i_binvar[binvar] = var_lims[binvar]

            output_dict = {"data_list":data_list_i_binvar, "file_list":[],"dir_list":[]}

            # Case that aggregate keys is length 0
            if len(aggregate_keys)==0:

                # Get job directory and output file name
                job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i[key])]) for key in sorted(data_list_i)]))
                job_dir = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
                out_file_name = get_out_file_name(use_mc=use_mc,binvar=binvar,binvar_min=binvar_min,binvar_max=binvar_max,**data_list_i)
                out_file_name = os.path.join(job_dir,out_file_name)
                output_dict["file_list"].append(out_file_name)
                output_dict["dir_list"].append(job_dir)
                
            # Loop aggregate keys and build file list for current binning
            for key in aggregate_keys:
                for value in divisions[key]:
                
                    data_list_i_val = dict(_data_list_i) #NOTE: CLONE OVERALL DICT NOT data_list_i SINCE THAT HAS BINNING VARIABLE LIMITS IN IT.
                    data_list_i_val[key] = value

                    # Get job directory and output file name
                    job_dir = os.path.join(base_dir,"__".join(["_".join([key,str(data_list_i_val[key])]) for key in sorted(data_list_i_val)]))
                    job_dir = os.path.abspath(job_dir) #NOTE: Since dictionary is not copied this should just edit the original entry in data_list.
                    out_file_name = get_out_file_name(binvar=binvar,binvar_min=binvar_min,binvar_max=binvar_max,**data_list_i)
                    out_file_name = os.path.join(job_dir,out_file_name)
                    output_dict["file_list"].append(out_file_name)
                    output_dict["dir_list"].append(job_dir)

            # Now add output_dict to your overall file list
            out_file_list.append(output_dict)

    return out_file_list

def offset_graph_x(g, offset):
    """
    Parameters
    ----------
    g : np.array, required
        Numpy array describing a graph with structure g[{x_mean,x_err,y_mean,y_err},point_idx]
    offset : float, required
        Value by which to offset graph x values

    Description
    -----------
    Offset the x values of a graph.
    """
    npoints = len(g[0])
    for idx in range(npoints):
        g[0][idx] += offset

def save_graph_to_csv(
        filename,
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
    if xerr is None or len(xerr)==0: xerr = [0.0 for el in x]
    if yerr is None or len(yerr)==0: yerr = [0.0 for el in x]
    if xerr_syst is None or len(xerr_syst)==0: xerr_syst = [0.0 for el in x]
    if yerr_syst is None or len(yerr_syst)==0: yerr_syst = [0.0 for el in x]
    for i, el in enumerate(x):
        data.append([i, x[i], y[i], xerr[i], yerr[i], xerr_syst[i], yerr_syst[i]])
    data = np.array(data)

    # Save data to file
    header = "REPLACEMENT_HEADER"+header
    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(filename, 'w') as file:
        file.write(filedata)

def save_systematics_to_csv(
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
    header = "REPLACEMENT_HEADER"+header
    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(filename, 'w') as file:
        file.write(filedata)

def load_TGraphError(
        path,
        name = "Graph"
    ):
    """
    Parameters
    ----------
    path : required, string
        Path to ROOT file containing TGraphError
    name : optional, string
        Name of TGraphError object with ROOT file

    Returns
    -------
    List
        Graph data as a list of numpy arrays with structure [`x`, `y`, `xerr`, `yerr`] or empty list if file is not found

    Description
    -----------
    Read TGraphError data from a ROOT file.
    """

    # Get TGraphErrors from ROOT file
    try:
        f = ur.open(path)
        g = [np.array(f[name].member('fX')), f[name].member('fY'), f[name].member('fEX'), f[name].member('fEY')]
        return g

    except FileNotFoundError as e:
        print("FileNotFoundError: ",path)
        print("Returning empty list")
        return []

def get_tgrapherror_arrs(
        out_file_list,
        name='Graph',
        sgasym=0.0
    ):
    """
    Parameters
    ----------
    out_file_list : required, list
        List of paths to output ROOT files containing graph data in TGraphError objects
    name : optional, string
        Name of TGraphError object with ROOT file
    sgasym : optional, float
        Injected signal asymmetry for computing difference of measured and injected values

    Returns
    -------
    Dictionary
        Map of mean x and y values and other statistics names to a list of their values in each kinematic bin

    Description
    -----------
    Compute the average x and y graph values as well as errors and other statistical information across a list of input files.
    The keys of the returned map are as follows:
    'x_mean'
    'y_mean'
    'xerr_mean'
    'yerr_mean'
    'y_min'
    'y_max'
    'y_std'
    'ydiff_mean'
    'ydiff_std'
    'ydiff_mins'
    'ydiff_maxs'
    """

    # Initialize output list
    glist = []
    
    # Loop files and load data from ROOT TGraphError objects
    for filename in out_file_list:
        g = load_TGraphError(filename,name=name)
        if len(g)>0: glist.append(g)

    # Check if graph list is empty
    if len(glist)==0:
        print("ERROR: len(glist)==0")
        return {
            'x_mean':[],
            'y_mean':[],
            'xerr_mean':[],
            'yerr_mean':[],
            'y_min':[],
            'y_max':[],
            'y_std':[],
            'ydiff_mean':[],
            'ydiff_std':[],
            'ydiff_mins':[],
            'ydiff_maxs':[],
            }

    # Convert to numpy and swap bin and graph axes
    glist  = np.array(glist)
    glist  = np.swapaxes(glist,0,1)

    # Get arrays
    x_mean     = np.mean(glist[0],axis=0) #NOTE: Get mean across different graphs (axis=0) but not across bins (axis=1)
    y_mean     = np.mean(glist[1],axis=0)
    xerr_mean = np.sqrt(np.mean(np.square(glist[2]),axis=0))
    yerr_mean = np.sqrt(np.mean(np.square(glist[3]),axis=0))
    y_min      = np.min(glist[1],axis=0)
    y_max      = np.max(glist[1],axis=0)
    y_std      = np.std(glist[1],axis=0)
    ydiff_mean = np.mean(glist[1]-sgasym,axis=0)
    ydiff_std  = np.std(glist[3]-sgasym,axis=0)
    ydiff_mins = np.min(glist[1]-sgasym,axis=0)
    ydiff_maxs = np.max(glist[1]-sgasym,axis=0)

    return {
            'x_mean':x_mean,
            'y_mean':y_mean,
            'xerr_mean':xerr_mean,
            'yerr_mean':yerr_mean,
            'y_min':y_min,
            'y_max':y_max,
            'y_std':y_std,
            'ydiff_mean':ydiff_mean,
            'ydiff_std':ydiff_std,
            'ydiff_mins':ydiff_mins,
            'ydiff_maxs':ydiff_maxs,
            }

def load_TH2(
        path,
        name = "h2"
    ):
    """
    Parameters
    ----------
    path : required, string
        Path to ROOT file containing TH2D
    name : optional, string
        Name of TH2D object with ROOT file

    Returns
    -------
    np.array
        Histogram data as a 2D numpy array or empty list if file is not found

    Description
    -----------
    Read TH2D data from a ROOT file.
    """

    # Get TH2D from ROOT file
    try:
        f = ur.open(path)
        g = f[name].values()
        return g

    except FileNotFoundError as e:
        print("FileNotFoundError: ",path)
        print("\t Returning empty list")
        return []

def save_bin_migration_matrix_to_csv(
        bin_migration_mat,
        base_dir='systematics/bin_migration/',
        binvar='Q2',
        delimiter=",",
        header=None,
        fmt=None,
        comments='',
    ):
    """
    Parameters
    ----------
    bin_migration_mat : required, np.array
        2D bin migration matrix with (i,j) `->` (generated,reconstructed)
    base_dir : string, required
        Path to directory in matrix will be saved
        Default : 'systematics/bin_migration/'
    binvar : optional, string
        Name of reconstructed bin variable
        Default : 'Q2'
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
    Save a 2D bin migration matrix mapping generated bins to reconstructed bins to a CSV file
    with an added initial row and column for the bin indices.  Note that files will be saved
    with name `bin_migration_mat_<binvar>.csv`.
    """

    if np.shape(bin_migration_mat)[0]!=np.shape(bin_migration_mat)[1] or len(np.shape(bin_migration_mat))!=2:
        raise TypeError("Bin migration matrix must be square but has shape "+str(np.shape(bin_migration_mat)))

    # Set output filename
    filename = 'bin_migration_mat_'+binvar+'.csv'
    filename = os.path.join(base_dir,filename)

    # Create new table with int bin labels
    nbins = np.shape(bin_migration_mat)[0]
    new_shape = list(np.shape(bin_migration_mat)) #NOTE: List is important here!
    new_shape[1] += 1 #NOTE: Increase the number of columns to accomodate bin numbers in the initial column
    data = np.zeros(new_shape)
    data[:,0] = [i for i in range(1,nbins+1)]
    data[0:,1:] = bin_migration_mat

    # Save to CSV
    if header is None: header = ' '+delimiter+delimiter.join([str(i) for i in range(1,nbins+1)])
    header = "REPLACEMENT_HEADER"+header
    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(filename, 'w') as file:
        file.write(filedata)

def load_csv(
        path,
        results_dir='results/',
        base_dir='systematics/'
    ):
    """
    Parameters
    ----------
    path : required, string
        Path to CSV file containing table
    results_dir : optional, string
        Directory name to replace
        Default : 'results/'
    base_dir : optional, string
        Directory name to insert
        Default : 'systematics/'

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe of csv file data

    Description
    -----------
    Read table stored in CSV format into a pandas DataFrame optionally replacing
    part of the path, for example, a directory name, with another value.
    """
    inpath = path.replace(results_dir,base_dir) if len(base_dir)>0 else path
    return pd.read_csv(inpath)

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

def plot_systematics(
        x_means,
        yerr_syst,
        palette = 'Dark2',
        stacked = False,
        label   = None,
        xlims   = (0.0,1.0),
        ylims   = (-1.0,1.0),
        xvar    = 'Q2',
        title   = 'Systematic Errors',
        xtitle  = '$Q^{2} (GeV^{2})$',
        ytitle  = '$\Delta \mathcal{A}$',
        outpath = 'systematics.pdf',
        add_clas12_watermark = True
    ):
    """
    Parameters
    ----------
    x_means : list, required
        Mean x values for each bin
    yerr_syst : np.array, required
        Array of absolute systematic error in each bin further indexed by the sources of systematic error
    palette : string, optional
        Default : 'Dark2'
    stacked : Boolean, optional
        Default : False
    label : List, optional
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
        Default : 'Q2'
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
    add_clas12_watermark : Boolean, optional
        Option to add CLAS12 watermark on produced plot
        Default : True

    Description
    -----------
    Plot the systematic error for each bin in a 1D binning scheme broken down by sources of systematic error.
    """

    # Set color palette
    sbn.set_palette(palette)

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

    # Set plotting parameters
    ecolor='black'
    elinewidth=2.0
    capsize=18
    capthick=2.0
    marker='o'
    linestyle=None
    linewidth=0.0
    markersize=20
    gridlinewidth=0.5
    axlinewidth=1

    # Set up plot
    figsize = (16,10)
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title,usetex=True)
    plt.xlabel(xtitle,usetex=True)
    plt.ylabel(ytitle,usetex=True)

    # Plot systematics by source for each x point
    xbins = x_means
    nbins = len(xbins)
    xbins = np.moveaxis(np.array([xbins for el in range(np.shape(yerr_syst)[1])]),(0,1),(1,0))
    s1 = plt.hist(xbins, weights=yerr_syst, bins=nbins, alpha=0.5, label=label, stacked=stacked)
    plt.tick_params(direction='out',bottom=True,top=True,left=True,right=True,length=10,width=1)
    ax1.axhline(0, color='black',linestyle='-',linewidth=axlinewidth)
    if add_clas12_watermark:
        plt.text(0.5, 0.5, 'CLAS12 Preliminary',
                size=50, rotation=25., color='gray', alpha=0.25,
                horizontalalignment='center',verticalalignment='center',transform=ax1.transAxes)
    plt.legend(loc='best')
    f1.savefig(outpath)
