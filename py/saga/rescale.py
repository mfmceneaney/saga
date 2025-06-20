"""
This module implements methods for rescaling asymmetry
results for future experiments.

# Author: Matthew F. McEneaney (2024, Duke University)
"""

import numpy as np
from .data import load_csv
from .aggregate import get_aggregate_graph


def rescale_graph_data(
    ct_mean,
    x_mean,
    y_mean,
    xerr_mean,
    yerr_mean,
    path,
    old_dat_path="old_dat_path.csv",
    new_sim_path="new_sim_path.csv",
    old_sim_path="old_sim_path.csv",
    count_key="count",
    yerr_key="",
    xs_ratio=1.0,
    lumi_ratio=1.0,
    tpol_factor=1.0,
    tdil_factor=1.0,
    yvalue=-100.0,
    xvar_keys=None,
    sgasym=0.0,
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
    new_sim_graph = load_csv(
        path, old_path=old_dat_path, new_path=new_sim_path, aliases=aliases
    )
    old_sim_graph = load_csv(
        path, old_path=old_dat_path, new_path=old_sim_path, aliases=aliases
    )

    # Get counts OR y errors from csv
    new_sim_graph_count = (
        new_sim_graph[count_key]
        if yerr_key is None or yerr_key == ""
        else 1.0 / np.square(new_sim_graph[yerr_key])
    )
    old_sim_graph_count = (
        old_sim_graph[count_key]
        if yerr_key is None or yerr_key == ""
        else 1.0 / np.square(old_sim_graph[yerr_key])
    )

    # Compute scaled quantities
    acceptanceratio = np.divide(new_sim_graph_count, old_sim_graph_count) / xs_ratio
    acceptanceratio[np.isinf(acceptanceratio)] = 0
    acceptanceratio[np.isnan(acceptanceratio)] = (
        0  # NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    )
    scaling = acceptanceratio * lumi_ratio
    scaled_ct_mean = np.multiply(scaling, ct_mean)
    err_scaling = np.sqrt(
        np.divide(ct_mean, scaled_ct_mean)
    )  # NOTE: SCALE ERRORS ASSUMING POISSONIAN STATISTICS -> d ~ 1/sqrt(N) -> MULTIPLY BY sqrt(N_old_data/N_new_data)
    err_scaling[np.isinf(err_scaling)] = 0
    err_scaling[np.isnan(err_scaling)] = (
        0  # NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    )
    scaled_yerr_mean = (
        np.multiply(err_scaling, yerr_mean) * 1.0 / (tpol_factor * tdil_factor)
    )

    # Set other graph quantities to new_sim values
    x_mean = new_sim_graph["x"]
    xerr_mean = new_sim_graph["xerr"]
    y_mean = new_sim_graph["y"]

    # Set y values to constant and update scaled y errors if requested
    scaled_y_mean = y_mean if yvalue < -1 else [yvalue for i in range(len(y_mean))]
    if yvalue >= -1:
        scaled_yerr_mean *= np.sqrt(1 - np.square(yvalue * tpol_factor))

    # Create a length 1 list of graph data with scaled graph results
    graph_list = np.array(
        [[scaled_ct_mean, scaled_y_mean, scaled_yerr_mean, x_mean, xerr_mean]]
    )

    graph = get_aggregate_graph(graph_list, xvar_keys=xvar_keys, sgasym=sgasym)

    # Set systematic errors to scaling fractions
    graph["scaling"] = scaling
    graph["acceptanceratio"] = acceptanceratio

    return graph


def rescale_csv_data(
    path,
    outpath="",
    old_dat_path="old_dat_path.csv",
    new_sim_path="new_sim_path.csv",
    old_sim_path="old_sim_path.csv",
    count_key="count",
    y_key="a0",
    yerr_key="a0err",
    xs_ratio=1.0,
    lumi_ratio=1.0,
    tpol_factor=1.0,
    tdil_factor=1.0,
    yvalue=-100.0,
    float_format="%.3g",
    config=None,
    aggregate_config=None,
    chain_configs=None,
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
    old_dat_df = load_csv(
        path,
        config=config,
        aggregate_config=aggregate_config,
        chain_configs=chain_configs,
        aliases=aliases,
    )
    new_sim_df = load_csv(
        path, old_path=old_dat_path, new_path=new_sim_path, aliases=aliases
    )
    old_sim_df = load_csv(
        path, old_path=old_dat_path, new_path=old_sim_path, aliases=aliases
    )

    # Get counts OR y errors from csv
    new_sim_df_count = new_sim_df[
        count_key
    ]  # if yerr_key is None or yerr_key == '' else 1.0/np.square(new_sim_df[yerr_key])
    old_sim_df_count = old_sim_df[
        count_key
    ]  # if yerr_key is None or yerr_key == '' else 1.0/np.square(old_sim_df[yerr_key])
    old_dat_df_count = old_dat_df[
        count_key
    ]  # if yerr_key is None or yerr_key == '' else 1.0/np.square(old_dat_df[yerr_key])

    # Compute scaled quantities
    yerrs = old_dat_df[yerr_key]
    acceptanceratio = np.divide(new_sim_df_count, old_sim_df_count) / xs_ratio
    acceptanceratio[np.isinf(acceptanceratio)] = 0
    acceptanceratio[np.isnan(acceptanceratio)] = (
        0  # NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    )
    scaling = acceptanceratio * lumi_ratio
    new_dat_df_count = np.multiply(scaling, old_dat_df_count)
    err_scaling = np.sqrt(
        np.divide(old_dat_df_count, new_dat_df_count)
    )  # NOTE: SCALE ERRORS ASSUMING POISSONIAN STATISTICS -> d ~ 1/sqrt(N) -> MULTIPLY BY sqrt(N_old_data/N_new_data)
    err_scaling[np.isinf(err_scaling)] = 0
    err_scaling[np.isnan(err_scaling)] = (
        0  # NOTE: Check for all cases of zero division (n>1/0 or n==0/0) and replace with zero
    )
    scaled_yerrs = np.multiply(err_scaling, yerrs) * 1.0 / (tpol_factor * tdil_factor)

    # Set y values to constant and update scaled y errors if requested
    scaled_ys = (
        old_dat_df[y_key]
        if yvalue < -1
        else [yvalue for i in range(len(old_dat_df[y_key]))]
    )
    if yvalue >= -1:
        scaled_yerrs *= np.sqrt(1 - np.square(yvalue * tpol_factor))

    # Copy the old dataframe into the new dataframe
    new_dat_df = old_dat_df.copy(deep=True)

    # Set the entries of your new dataframe
    new_dat_df[count_key] = new_dat_df_count
    new_dat_df[y_key] = scaled_ys
    new_dat_df[yerr_key] = scaled_yerrs
    new_dat_df["scaling"] = scaling
    new_dat_df["acceptanceratio"] = acceptanceratio

    # Save the new dataframe to CSV
    outpath = path.replace(".csv", "_rescaled.csv") if len(outpath) == 0 else outpath
    new_dat_df.to_csv(outpath, float_format=float_format, index=False)
