"""
This module implements methods for loading and saving outputs
from slurm jobs for all possible option combinations
supplied in yaml files.

# Author: Matthew F. McEneaney (2024, Duke University)
"""

import os
import yaml
import uproot as ur
import numpy as np
import pandas as pd
from .aggregate import get_config_list, get_config_str


def load_th1(path, name="h1"):
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
    Read :obj:`TH1` histogram data from a ROOT file.
    This will work for any histogram: (:obj:`TH1`, :obj:`TH2`, :obj:`TH3`).
    """

    # Get TH1 from ROOT file
    if path == "":
        return []
    try:
        f = ur.open(path)
        g = f[name].to_numpy()
        return g

    except FileNotFoundError:
        print("FileNotFoundError: ", path)
        print("\t Returning empty list")
        return []


def load_yaml(path):
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
    with open(path, encoding="utf-8") as f:
        yaml_args = yaml.safe_load(f)
    return yaml_args


def load_csv(
    path,
    old_path=None,
    new_path=None,
    config=None,
    aggregate_config=None,
    chain_configs=None,
    aliases=None,
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
        Map of aggregate configuration option names to option values,
        used for determining correct directory names
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
    # Set input path
    inpath = path
    if old_path is not None and new_path is not None:
        inpath = path.replace(old_path, new_path)

    # Chain CSVs across the given keys
    if (
        config is not None
        and len(config) > 0
        and chain_configs is not None
        and len(chain_configs) > 0
    ):

        # Get the full batch config
        aggregate_configs = {}
        if aggregate_config is not None:
            aggregate_configs = {
                key: [aggregate_config[key]] for key in aggregate_config
            }
        configs = dict(
            {key: [config[key]] for key in config}, **chain_configs, **aggregate_configs
        )

        # Get a list of all possible option value combinations from configs
        config_list = get_config_list(configs, aggregate_keys=[])

        # Set csv list
        csv_list = []

        # Loop resulting list
        for config_list_i in config_list:

            # Get job directory
            config_str = get_config_str(config_list_i, aliases=aliases)

            # Get base job directory
            base_config_str = get_config_str(config, aliases=aliases)

            # Modify path for chain element
            inpath_i = inpath.replace(base_config_str, config_str)

            # Open csv
            csv_i = pd.read_csv(inpath_i)
            csv_list.append(csv_i)

        # Merge csvs and return
        return pd.concat(csv_list, ignore_index=True)

    # Return csv
    return pd.read_csv(inpath)


def save_txt(
    filename,
    data,
    delimiter=",",
    header=None,
    fmt=None,
    comments="",
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
    if header is None:  # NOTE: ASSUME DATA HAS DIMENSION: [NCOL,NROWS]
        header = delimiter.join([str(i) for i in range(len(data))])
    if fmt is None:
        fmt = delimiter.join(["%.3g" for el in data])
    np.savetxt(
        filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments
    )


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
    comments="",
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
        Raise an error if the shape of the systematics arrays
        do not match :obj:`(nbins)` or :obj:`(nbins,2)`.


    Description
    -----------
    Write a graph to a CSV file with optional errors and systematic errors.
    Systematic errors may have high and low values.
    """

    # Create data array
    data = []
    if ct is None or len(ct) == 0:
        ct = [0.0 for el in x]
    if xerr is None or len(xerr) == 0:
        xerr = [0.0 for el in x]
    if yerr is None or len(yerr) == 0:
        yerr = [0.0 for el in x]
    if xerr_syst is None or len(xerr_syst) == 0:
        xerr_syst = [0.0 for el in x]
    if yerr_syst is None or len(yerr_syst) == 0:
        yerr_syst = [0.0 for el in x]
    xerr_syst_shape = np.shape(xerr_syst)
    yerr_syst_shape = np.shape(yerr_syst)
    for i, el in enumerate(x):
        data_i = [i, ct[i], x[i], y[i], xerr[i], yerr[i]]
        if len(xerr_syst_shape) == 1 and len(yerr_syst_shape) == 1:
            data_i.extend([xerr_syst[i], yerr_syst[i]])
        elif (
            len(xerr_syst_shape) == 2
            and xerr_syst_shape[1] == 2
            and len(yerr_syst_shape) == 1
        ):
            data_i.extend([xerr_syst[i][0], xerr_syst[i][1], yerr_syst[i]])
        elif (
            len(xerr_syst_shape) == 1
            and len(yerr_syst_shape) == 2
            and yerr_syst_shape[1] == 2
        ):
            data_i.extend([xerr_syst[i], yerr_syst[i][0], yerr_syst[i][1]])
        elif (
            len(xerr_syst_shape) == 2
            and xerr_syst_shape[1] == 2
            and len(yerr_syst_shape) == 2
            and yerr_syst_shape[1] == 2
        ):
            data_i.extend(
                [xerr_syst[i][0], xerr_syst[i][1], yerr_syst[i][0], yerr_syst[i][1]]
            )
        else:
            raise ValueError(
                f"ERROR: xerr_syst_shape={xerr_syst_shape} or yerr_syst_shape={yerr_syst_shape} "
                + "does not have shape (nbins) or (nbins,2)."
            )
        data.append(data_i)
    data = np.array(data)

    # Save data to file
    save_txt(
        filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments
    )


def save_graph_systematics_to_csv(
    filename, x, yerrs_syst=None, delimiter=",", header=None, fmt=None, comments=""
):
    """
    Parameters
    ----------
    filename : str, required
        Output file name
    x : list, required
        Graph x values with shape :obj:`(nbins)`
    yerr_syst : list, optional
        Graph y systematic error values decomposed into the different sources of
        systematic error with shape :obj:`(nbins,nsources)` or :obj:`(nbins,nsources,2)`
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
        Raise an error if the shape of the systematics is does not match
        :obj:`(nbins,nsources)` or :obj:`(nbins, nsources,2)`.

    Description
    -----------
    Write a set of graph y systematic errors to a CSV file with the
    systematic error values broken down by source and allowing high and low errors.
    This means the argument :obj:`yerr_syst` should have shape :obj:`(nbins, nsources)`
    or :obj:`(nbins, nsources,2)` where :obj:`nbins` is the number of kinematic bins
    and :obj:`nsources` is the number of sources of systematic error.
    """

    # Create data array
    data = []
    if yerrs_syst is None or len(yerrs_syst) == 0:
        yerrs_syst = [[0.0] for el in x]
    yerrs_syst_shape = np.shape(yerrs_syst)
    if yerrs_syst_shape[0] == len(x) and len(yerrs_syst_shape) == 2:
        for i, el in enumerate(x):
            data.append([i, el, *yerrs_syst[i]])
    elif yerrs_syst_shape[0] == len(x) and (
        len(yerrs_syst_shape) == 3 and yerrs_syst_shape[2] == 2
    ):
        for i, el in enumerate(x):
            data_i = [i, el]
            for source in yerrs_syst[i]:
                data_i.extend(*source)
            data.append(data_i)
    else:
        raise ValueError(
            f"ERROR: yerrs_syst has shape {yerrs_syst} "
            + f"but allowed shapes are ({len(x)},*) and ({len(x)},*,2)."
        )
    data = np.array(data)

    # Save data to file
    save_txt(
        filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments
    )


def save_bin_mig_mat_to_csv(
    bin_mig_mat,
    base_dir="./",
    basename="",
    delimiter=",",
    header=None,
    fmt=None,
    comments="",
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
    if (
        np.shape(bin_mig_mat)[0] != np.shape(bin_mig_mat)[1]
        or len(np.shape(bin_mig_mat)) != 2
    ):
        raise TypeError(
            "Bin migration matrix must be square but has shape "
            + str(np.shape(bin_mig_mat))
        )

    # Set output filename
    filename = "bin_mig_mat_" + basename + ".csv"
    filename = os.path.join(base_dir, filename)

    # Create new table with int bin labels
    nbins = np.shape(bin_mig_mat)[0]
    new_shape = list(np.shape(bin_mig_mat))  # NOTE: List is important here!
    new_shape[1] += 1  # NOTE: Add bin numbers.
    data = np.zeros(new_shape)
    data[:, 0] = list(range(1, nbins + 1))
    data[0:, 1:] = bin_mig_mat

    # Set column formats if not given
    if fmt is None:
        fmt = ["%.3g" for i in range(np.shape(bin_mig_mat)[0])]
        fmt = ["%d", *fmt]

    save_txt(
        filename, data, header=header, delimiter=delimiter, fmt=fmt, comments=comments
    )
