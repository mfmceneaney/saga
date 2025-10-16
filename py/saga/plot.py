"""
This module implements methods for plotting
asymmetry results.
"""

import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import seaborn as sbn
from .data import load_th1, save_graph_to_csv, save_graph_systematics_to_csv
from .rescale import rescale_graph_data


def set_default_plt_settings():
    """
    Description
    -----------
    Set plt.rc parameters for font sizes and family and tick font size and tick length and direction
    in a nice format.
    """

    # Use LaTeX for text rendering
    plt.rcParams['text.usetex'] = True

    # Set font sizes
    plt.rc("font", size=25)  # controls default text size
    plt.rc("axes", titlesize=50)  # fontsize of the title
    plt.rc("axes", labelsize=50)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=25)  # fontsize of the x tick labels
    plt.rc("ytick", labelsize=25)  # fontsize of the y tick labels
    plt.rc("legend", fontsize=25)  # fontsize of the legend

    # Get some nicer plot settings
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["figure.autolayout"] = True

    # Set tick parameters
    plt.tick_params(
        direction="out",
        bottom=True,
        top=True,
        left=True,
        right=True,
        length=10,
        width=1,
    )


def plot_injected_asyms(
    ax1,
    asyms,
    ytitles,
    colors,
    sgasym_idx=0,
    ylims=(-1.0, 1.0),
    label_base="Injected Signal ",
    linestyle="--",
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
    plotted_values = (
        {}
    )  # NOTE: Keep track of how many times you've plotted each constant asymmetry value
    for idx, asym_i in enumerate(asyms):

        # Plot injected asymmetries as (x,y) data OR axis lines
        if len(np.shape(asyms)) > 1:
            ax1.plot(
                asym_i[0],
                asym_i[1],
                color=colors[idx] if idx < len(colors) else None,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=0.5 if idx != sgasym_idx else 1.0,
                label=label_base + ytitles[idx] if idx < len(ytitles) else None,
            )
        else:
            # Set the offset
            offset = 0.0
            if asym_i in plotted_values:
                offset = plotted_values[asym_i] * 0.0025 * (ylims[1] - ylims[0])
                plotted_values[asym_i] += 1
            else:
                plotted_values[asym_i] = 1

            # Flip offset if the asymmetry is negative
            if asym_i < 0.0:
                offset *= -1.0

            # Plot the constant asymmetry
            ax1.axhline(
                asym_i + offset,
                color=colors[idx] if idx < len(colors) else None,
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=0.5 if idx != sgasym_idx else 1.0,
                label=label_base + ytitles[idx] if idx < len(ytitles) else None,
            )


def plot_watermark(ax1, watermark="CLAS12 Preliminary"):
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
    plt.text(
        0.5,
        0.5,
        watermark,
        size=50,
        rotation=25.0,
        color="gray",
        alpha=0.25,
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax1.transAxes,
    )


def plot_vlines(
    hist,
    binlims=None,
    linestyle="dotted",
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

    # Check arguments
    if binlims is None:
        binlims = []

    # Loop middle bin limits
    for xval in binlims[1:-1]:

        # Loop histogram x values
        for idx in range(len(hist[1]) - 1):

            # Check if bin limit is in bin
            binx_low, binx_high = hist[1][idx], hist[1][idx + 1]
            if binx_low <= xval < binx_high:

                # Plot lower bin limit
                ymax = hist[0][idx]
                plt.vlines(xval, 0.0, ymax, linestyle=linestyle)


def plot_hists(
    ax1,
    hist_paths=None,
    hist_keys=None,
    clone_axis=True,
    ylabel="Density",
    ylims=(0.0, 0.05),
    histtype="step",
    hist_colors=None,
    alpha=0.5,
    linewidth=2,
    density=True,
    log=False,
    hist_labels=None,
    binlims=None,
    vlinestyle="dotted",
    vline_hist_idx=-1,
    legend_loc="upper right",
    hist_dim=1,
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

    # Check arguments
    if hist_paths is None:
        hist_paths = []
    if hist_keys is None:
        hist_keys = []
    if binlims is None:
        binlims = []

    # Clone y axis and set labels
    ax2 = (
        ax1.twinx() if clone_axis else ax1
    )  # instantiate a second y-axis that shares the same x-axis
    if clone_axis:
        ax2.set_ylabel(ylabel)
        ax2.set_ylim(*ylims)

    # Set colors to just be taken from the current color palette if not provided
    if hist_colors is None:
        hist_colors = [None for i in range(len(hist_paths))]

    # Set histogram labels to just be taken from the current color palette if not provided
    if hist_labels is None:
        hist_labels = ["h" + str(i) for i in range(len(hist_paths))]

    # Check line width as a flag for drawing th2ds
    if hist_dim == 1:

        # Loop histograms and draw
        for idx, hist_path in enumerate(hist_paths):

            # Load histogram and convert to numpy
            h_y, h_bins = load_th1(hist_path, hist_keys[idx])

            # Get mean x bin values
            h_x = [(h_bins[i] + h_bins[i + 1]) / 2 for i in range(len(h_bins) - 1)]

            # Plot histogram
            h = ax2.hist(
                h_x,
                bins=h_bins,
                weights=h_y / np.sum(h_y) if density else h_y,
                histtype=histtype,
                color=hist_colors[idx],
                alpha=alpha,
                linewidth=linewidth,
                label=hist_labels[idx],
                density=False,
                log=log,
            )

            # Plot bin limits if supplied and we are on the last histogram
            if (
                idx == (vline_hist_idx if vline_hist_idx >= 0 else len(hist_paths) - 1)
                and len(binlims) > 0
            ):
                plot_vlines(
                    h,
                    binlims,
                    linestyle=vlinestyle,
                )

    # Assume TH2D histograms and plot
    else:

        # Loop histograms and draw
        for idx, hist_path in enumerate(hist_paths):

            # Load histogram and convert to numpy
            h2 = load_th1(hist_path, hist_keys[idx])

            # Plot histogram
            plot_th2(h2, ax1, add_colorbar=True, norm=LogNorm(), label=hist_labels[idx])

    # Plot legend if you cloned axis
    if clone_axis and legend_loc is not None and legend_loc != "":
        ax2.legend(loc=legend_loc)


def get_bin_kinematics_title(
    bin_id, df, cols=None, col_titles=None, err_ext="_err", sep=" , "
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

    # Check arguments
    if cols is None:
        cols = []
    if col_titles is None:
        col_titles = {}

    return sep.join(
        [
            f"$<{col_titles[col]}> = {df.iloc[bin_id].loc[col]:.2f}\\pm{df.iloc[bin_id].loc[col+err_ext]:.2f}$"
            for idx, col in enumerate(cols)
        ]
    )


def get_lims_coords(
    node,
    outer_xlims,
    outer_ylims,
    var_keys=None,
    nested_key="nested",
    lims_key="lims",
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

    # Initialize coordinates list and variable keys list
    lims_coords = []
    if var_keys is None:
        var_keys = []

    # Check node
    if (
        isinstance(node, dict)
        and nested_key in node
        and isinstance(node[nested_key], list)
        and isinstance(node[nested_key][0], dict)
    ):

        # Get first level nested node and get horizontal limits coordinates
        binvar_x = list(node[nested_key][0].keys())[0]
        node_nested = node[nested_key][0][binvar_x]
        horizontal_lims = [
            [outer_xlims, [y0, y0]] for y0 in node_nested[lims_key][1:-1]
        ]
        ylims = node_nested[lims_key]

        # Check nested node and get limits list
        xlims = []
        if (
            isinstance(node_nested, dict)
            and nested_key in node_nested
            and isinstance(node_nested[nested_key], list)
            and isinstance(node_nested[nested_key][0], dict)
        ):

            # Loop nested nodes and get limits lists
            for el in node_nested[nested_key]:
                binvar_y = list(el.keys())[0]
                if lims_key in el[binvar_y]:
                    xlims.append(el[binvar_y][lims_key])

            # Loop xlims and set vertical limit coordinates
            vertical_lims = []
            for yidx, xlim in enumerate(xlims):
                for xidx in range(len(xlim[1:-1])):
                    el = [
                        [xlim[xidx + 1], xlim[xidx + 1]],
                        [
                            ylims[yidx] if yidx > 0 else outer_ylims[0],
                            (
                                ylims[yidx + 1]
                                if yidx < len(xlims) - 1
                                else outer_ylims[-1]
                            ),
                        ],
                    ]
                    vertical_lims.append(el)

            # Add to coordinates list
            lims_coords.extend(vertical_lims)
        lims_coords.extend(horizontal_lims)

    # Grid scheme case
    elif isinstance(node, dict) and len(var_keys) == 2:

        # Get limits from node
        xvar, yvar = var_keys
        xlims = node[xvar]
        ylims = node[yvar]

        # Get line coordinates
        horizontal_lims = [[outer_xlims, [y0, y0]] for y0 in ylims[1:-1]]
        vertical_lims = [[[x0, x0], outer_ylims] for x0 in xlims[1:-1]]

        # Add to coordinates list
        lims_coords.extend(vertical_lims)
        lims_coords.extend(horizontal_lims)

    # Default to raising value error
    else:
        raise ValueError(
            f"Could not identify a 2D nested or grid bin scheme in `node`:\n{node}"
        )

    # Swap axes if needed
    if swap_axes:
        lims_coords = [[el[1], el[0]] for el in lims_coords]

    return lims_coords


def plot_lines(ax, coordinates, linecolor="red", linewidth=1):
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
    if len(np.shape(coordinates)) == 3 and np.shape(coordinates)[1:] != (2, 2):
        raise ValueError(f"Expected shape (nlines,2,2) but got {np.shape(coordinates)}")

    # Plot lines
    for coords in coordinates:
        ax.plot(*coords, color=linecolor, linewidth=linewidth, marker="o", markersize=0)


def get_bin_centers(cuts, swap_axes=False):
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
    new_cuts = list(cuts.values())
    bin_ids = list(cuts.keys())

    # Convert 2D bin cuts into bin widths and centers
    bin_centers = [cut.replace("=", "").split(") && (") for cut in new_cuts]
    bin_centers = [
        [el[0].replace("(", "").split(" && "), el[1].replace(")", "").split(" && ")]
        for el in bin_centers
    ]
    bin_centers = [
        [
            [float(el_idx[0].split(">")[1]), float(el_idx[1].split("<")[1])]
            for idx, el_idx in enumerate(el)
        ]
        for el in bin_centers
    ]
    bin_widths = [
        [el[0][1] - el[0][0], el[1][1] - el[1][0]] for el in bin_centers
    ]  # NOTE: ORDERING IS IMPORTANT HERE
    bin_centers = [
        [np.average(el[0]).item(), np.average(el[1]).item()] for el in bin_centers
    ]

    # Swap axes if needed
    if swap_axes:
        bin_centers = [[el[1], el[0]] for el in bin_centers]
        bin_widths = [[el[1], el[0]] for el in bin_widths]

    # Convert back to dictionaries
    bin_centers = {bin_ids[idx]: bin_centers[idx] for idx in range(len(bin_centers))}
    bin_widths = {bin_ids[idx]: bin_widths[idx] for idx in range(len(bin_widths))}

    return bin_centers, bin_widths


def plot_bin_ids(
    ax,
    bin_centers,
    size=25,
    color="red",
    alpha=1.0,
):
    """
    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes, required
        Matplotlib.pyplot axis to plot on
    bin_centers : dict, required
        Dictionary of bin ids to bin centers
    size : int, optional
        Font size for bin id text
    color : str, optional
        Text color
    alpha : float, optional
        Text alpha value
    """

    # Plot bin ids on bin centers
    for bin_id in bin_centers:
        ax.text(
            *bin_centers[bin_id],
            f"{bin_id}",
            size=size,
            color=color,
            alpha=alpha,
            horizontalalignment="center",
            verticalalignment="center",
        )


def plot_th2(h2, ax, add_colorbar=True, norm=LogNorm(), **kwargs):
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
    x = np.ravel(
        [
            [np.average([h2[1][i], h2[1][i + 1]]) for j in range(len(h2[2]) - 1)]
            for i in range(len(h2[1]) - 1)
        ]
    )
    y = np.ravel(
        [
            [np.average([h2[2][j], h2[2][j + 1]]) for j in range(len(h2[2]) - 1)]
            for i in range(len(h2[1]) - 1)
        ]
    )

    # Get the counts in each bin
    weights = np.ravel(h2[0])

    # Get the bin sizes
    bins = (h2[1], h2[2])

    # Plot the histogram
    hist2d = ax.hist2d(x, y, weights=weights, bins=bins, norm=norm, **kwargs)
    if add_colorbar:
        plt.colorbar(hist2d[3], ax=ax)


def plot_systematics(
    x_means,
    yerr_syst,
    palette="Dark2",
    stacked=False,
    syst_names=None,
    syst_labels=None,
    xlims=(0.0, 1.0),
    ylims=(-1.0, 1.0),
    title="Systematic Errors",
    xtitle="$Q^{2} (GeV^{2})$",
    ytitle="$\\Delta \\mathcal{A}$",
    outpath="systematics.pdf",
    watermark="CLAS12 Preliminary",
    use_default_plt_settings=True,
    legend_loc="upper left",
    axlinewidth=1.0,
    log=False,
    figsize=(16, 10),
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
    if use_default_plt_settings:
        set_default_plt_settings()

    # Set up plot
    f1, ax1 = plt.subplots(figsize=figsize)
    plt.xlim(*xlims)
    plt.ylim(*ylims)
    plt.title(title, usetex=True)
    plt.xlabel(xtitle, usetex=True)
    plt.ylabel(ytitle, usetex=True)

    # Plot systematics by source for each x point
    nbins = len(x_means)
    xbins = np.moveaxis(
        np.array([x_means for el in range(np.shape(yerr_syst)[1])]), (0, 1), (1, 0)
    )
    plt.hist(
        xbins,
        weights=yerr_syst,
        bins=nbins,
        alpha=0.5,
        label=syst_labels,
        stacked=stacked,
        log=log,
    )

    # Plot zero line
    ax1.axhline(0, color="black", linestyle="-", linewidth=axlinewidth)

    # Add water mark
    if watermark is not None and watermark != "":
        plot_watermark(ax1, watermark=watermark)

    # Plot legend
    if legend_loc is not None and legend_loc != "":
        ax1.legend(loc=legend_loc)

    # Save figure
    f1.savefig(outpath)

    # Save plot data to csv
    delimiter = ","
    if syst_names is None:
        syst_names = ["syst" + str(idx) for idx in range(nbins)]
    header = delimiter.join(
        ["bin", *syst_names]
    )  # NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    syst_fmts = ["%.3g" for idx in range(nbins)]
    fmt = ["%d", *syst_fmts]
    comments = ""

    # Save to CSV
    save_graph_systematics_to_csv(
        outpath + ".csv",
        x_means,
        yerrs_syst=yerr_syst,
        delimiter=delimiter,
        header=header,
        fmt=fmt,
        comments=comments,
    )


def plot_results(
    ax1,
    ct_mean=None,
    x_mean=None,
    y_mean=None,
    xerr_mean=None,
    yerr_mean=None,
    xerr_syst=None,
    yerr_syst=None,
    y_std=None,
    ydiff_mean=None,
    ydiff_std=None,
    xlims=(0.0, 1.0),
    ylims=(-1.0, 1.0),
    title="Asymmetry Results",
    xvar="x",
    xlabel="$x$",
    ylabel="$\\mathcal{A}$",
    sgasym_labels=None,
    bgasym_labels=None,
    sgasym_idx=0,
    sgasyms=None,
    bgasyms=None,
    sg_colors=None,
    bg_colors=None,
    fill_color="gray",
    outpath="out.pdf",
    watermark="CLAS12 Preliminary",
    show_injected_asymmetries=False,
    legend_loc="upper left",
    ecolor="black",
    elinewidth=2.0,
    capsize=18,
    marker="o",
    markersize=20,
    linestyle=None,
    linewidth=0.0,
    axlinewidth=1.0,
    hist_paths=None,
    hist_keys=None,
    hist_clone_axis=True,
    hist_ylabel="Density",
    hist_ylims=(0.0, 0.05),
    histtype="step",
    hist_colors=None,
    hist_alpha=0.5,
    hist_linewidth=2,
    hist_density=True,
    hist_log=False,
    hist_labels=None,
    binlims=None,
    vlinestyle="dotted",
    vline_hist_idx=-1,
    hist_legend_loc="upper right",
    hist_dim=1,
    old_dat_path="old_dat_path.csv",
    new_sim_path="new_sim_path.csv",
    old_sim_path="old_sim_path.csv",
    count_key="count",
    yerr_key="",
    xs_ratio=1.0,
    lumi_ratio=0.0,
    tpol_factor=1.0,
    tdil_factor=1.0,
    graph_yvalue=-100.0,
    aliases=None,
    plot_xerrors=False,
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
    y_std : list, optional
        y standard deviation values for each bin with shape :obj:`(nbins)`
    ydiff_mean : list, optional
        y difference from injected signal asymmetry mean values for each bin with shape :obj:`(nbins)`
    ydiff_std : list, optional
        y difference from injected signal asymmetry standard deviation values for each bin with shape :obj:`(nbins)`
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
    marker : str, optional
        Marker type
    markersize : int, optional
        Marker size
    linestyle : str, optional
        Line style
    linewidth : float, optional
        Line width
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

    # Check arguments
    if sgasym_labels is None:
        sgasym_labels = ["$\\mathcal{A}$"]
    if sgasym_idx is None:
        sgasym_idx = 0
    if sgasyms is None:
        sgasyms = [0.10]
    if sg_colors is None:
        sg_colors = ["blue"]

    # Rescale graph
    scaling, acceptanceratio = None, None
    rescale = lumi_ratio > 0.0
    if rescale:

        # Get rescaled graph data
        rescaled_graph = rescale_graph_data(
            ct_mean,
            x_mean,
            y_mean,
            xerr_mean,
            yerr_mean,
            outpath + ".csv",
            old_dat_path=old_dat_path,
            new_sim_path=new_sim_path,
            old_sim_path=old_sim_path,
            count_key=count_key,
            yerr_key=yerr_key,
            xs_ratio=xs_ratio,
            lumi_ratio=lumi_ratio,
            tpol_factor=tpol_factor,
            tdil_factor=tdil_factor,
            yvalue=graph_yvalue,
            xvar_keys=[xvar],
            sgasym=sgasyms[sgasym_idx],
            aliases=aliases,
        )

        # Reset graph data
        ct_mean = rescaled_graph["ct_mean"]
        x_mean = rescaled_graph["x_mean"]
        y_mean = rescaled_graph["y_mean"]
        xerr_mean = rescaled_graph["xerr_mean"]
        yerr_mean = rescaled_graph["yerr_mean"]
        y_std = rescaled_graph["y_std"]
        ydiff_mean = rescaled_graph["ydiff_mean"]
        ydiff_std = rescaled_graph["ydiff_std"]
        scaling = rescaled_graph["scaling"]
        acceptanceratio = rescaled_graph["acceptanceratio"]

    # Set up plot
    ax1.set_xlim(*xlims)
    ax1.set_ylim(*ylims)
    ax1.set_title(title, usetex=True)
    ax1.set_xlabel(xlabel, usetex=True)
    ax1.set_ylabel(ylabel, usetex=True)

    # Plot projection variable distribution histograms
    if len(hist_paths) > 0:
        plot_hists(
            ax1,
            hist_paths=hist_paths,
            hist_keys=hist_keys,
            clone_axis=hist_clone_axis,
            ylabel=hist_ylabel,
            ylims=hist_ylims,
            histtype=histtype,
            hist_colors=hist_colors,
            alpha=hist_alpha,
            linewidth=hist_linewidth,
            density=hist_density,
            log=hist_log,
            hist_labels=hist_labels,
            binlims=binlims,
            vlinestyle=vlinestyle,
            vline_hist_idx=vline_hist_idx,
            legend_loc=hist_legend_loc,
            hist_dim=hist_dim,
        )

    # Plot systematic errors
    if yerr_syst is not None:
        ax1.errorbar(
            x_mean,
            y_mean,
            xerr=None,
            yerr=yerr_syst,
            ecolor=fill_color,
            elinewidth=elinewidth * 20,
            capsize=0,
            color=fill_color,
            marker="o",
            linestyle=linestyle,
            alpha=0.5,
            linewidth=0,
            markersize=0,
            label="Systematic error",
        )

    # Plot standard deviation of aggregated injected values
    if y_std is not None:
        ax1.fill_between(
            x_mean,
            np.add(y_mean, y_std),
            np.add(y_mean, -y_std),
            alpha=0.2,
            label="$\\pm1\\sigma$ Band",
            color=fill_color,
        )

    # Plot results
    if yerr_mean is not None:
        ax1.errorbar(
            x_mean,
            y_mean,
            xerr=xerr_mean if plot_xerrors else None,
            yerr=yerr_mean,
            ecolor=ecolor,
            elinewidth=elinewidth,
            capsize=capsize,
            color=sg_colors[sgasym_idx],
            marker=marker,
            linestyle=linestyle,
            linewidth=linewidth,
            markersize=markersize,
            label=sgasym_labels[sgasym_idx],
        )

    # Add zero line
    ax1.axhline(0, color="black", linestyle="-", linewidth=axlinewidth)

    # Draw injected asymmetries
    if show_injected_asymmetries:

        # Plot injected signal asymmetries
        plot_injected_asyms(
            ax1,
            sgasyms,
            sgasym_labels,
            sg_colors,
            sgasym_idx=sgasym_idx,
            ylims=ylims,
            label_base="Injected Signal ",
            linestyle="--",
            linewidth=axlinewidth,
        )

        # Plot injected background asymmetries
        plot_injected_asyms(
            ax1,
            bgasyms,
            bgasym_labels,
            bg_colors,
            sgasym_idx=-1,  # NOTE: Don't emphasize any background asymmetries for now.
            ylims=ylims,
            label_base="Injected Background ",
            linestyle="--",
            linewidth=axlinewidth,
        )

    # Add water mark
    if watermark is not None and watermark != "":
        plot_watermark(ax1, watermark=watermark)

    # Plot legend
    if legend_loc is not None and legend_loc != "":
        ax1.legend(loc=legend_loc)

    # Check whether you have graph data to save to CSV
    if ct_mean is None:
        return

    # Save plot data to csv
    delimiter = ","
    cols = (
        ["bin", "count", "x", "y", "xerr", "yerr", "xerrsyst", "yerrsyst"]
        if not rescale
        else ["bin", "count", "x", "y", "xerr", "yerr", "acceptanceratio", "scaling"]
    )
    xerr_syst_shape = np.shape(xerr_syst)
    yerr_syst_shape = np.shape(yerr_syst)
    if len(xerr_syst_shape) == 2:
        cols = [
            "bin",
            "count",
            "x",
            "y",
            "xerr",
            "yerr",
            "xerrsystlow",
            "xerrsysthigh",
            "yerrsyst",
        ]
    if len(yerr_syst_shape) == 2:
        cols = [
            "bin",
            "count",
            "x",
            "y",
            "xerr",
            "yerr",
            "xerrsyst",
            "xerrsyst",
            "yerrsystlow",
            "yerrsysthigh",
        ]
    if len(xerr_syst_shape) == 2 and len(yerr_syst_shape) == 2:
        cols = [
            "bin",
            "count",
            "x",
            "y",
            "xerr",
            "yerr",
            "xerrsystlow",
            "xerrsysthigh",
            "yerrsystlow",
            "yerrsysthigh",
        ]
    header = delimiter.join(
        cols
    )  # NOTE: CAN'T HAVE UNDERSCORE IN COLUMN NAMES FOR LATEX CSVSIMPLE
    fmt = ["%.3g" for i in range(len(cols) - 2)]
    fmt = ["%d", "%d", *fmt]
    comments = ""

    # Save plot data
    save_graph_to_csv(
        outpath + ".csv" if not rescale else outpath + "_rescaled.csv",
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
        comments=comments,
    )

    # Save ydiffs for MC asym injection systematics
    if ydiff_mean is not None:
        save_graph_to_csv(
            outpath + "_ydiff.csv" if not rescale else outpath + "_rescaled_ydiff.csv",
            ct_mean,
            x_mean,
            ydiff_mean,
            xerr=xerr_mean,
            yerr=ydiff_std,
            delimiter=delimiter,
            header=header,
            fmt=fmt,
            comments=comments,
        )


def plot_results_array(
    graph_array,
    plot_results_kwargs_array,
    plot_results_kwargs_base=None,
    figsize=(16, 10),
    outpath="plot_projections.pdf",
    use_grid_titles=True,
    use_grid_xlabels=True,
    use_grid_ylabels=True,
    use_grid_hist_ylabels=True,
    use_default_plt_settings=True,
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

    # Check arguments
    if plot_results_kwargs_base is None:
        plot_results_kwargs_base = {}

    # Use default plotting settings
    if use_default_plt_settings:
        set_default_plt_settings()

    # Get and check graph matrix shape
    shape = np.shape(graph_array)
    if len(shape) not in (1, 2):
        raise TypeError(
            "`plot_projections()` : `graph_array` shape must have shape with len(shape) in (1,2) but shape = ",
            shape,
        )

    # Create figure and axes
    f, ax = plt.subplots(
        *shape, figsize=figsize, squeeze=not len(shape) > 1
    )  # NOTE: squeeze will squeeze out dimension one axes!

    # Loop axes and plot results for 1D and 2D cases
    if len(shape) == 1:
        for i in range(shape[0]):

            # Check for masked entry
            if graph_array[i] is None:
                continue

            # Format graph titles and axes depending on location in grid array
            if i != 0 and use_grid_titles:
                plot_results_kwargs_array[i]["title"] = ""
            if i != shape[0] - 1 and use_grid_xlabels:
                plot_results_kwargs_array[i]["xlabel"] = ""

            # Plot results
            plot_results_kwargs = dict(
                plot_results_kwargs_base, **plot_results_kwargs_array[i]
            )
            outpath_i = outpath.split(".")
            ext = outpath_i.pop(-1)
            outpath_i = ".".join(outpath_i)
            outpath_i = "".join([outpath_i, f"___arrloc_{i}.", ext])
            plot_results_kwargs["outpath"] = outpath_i
            plot_results(ax[i], **graph_array[i], **plot_results_kwargs)
    else:
        for i in range(shape[0]):
            for j in range(shape[1]):

                # Check for masked entry
                if graph_array[i][j] is None:
                    continue

                # Format graph titles and axes depending on location in grid array
                if j != 0 and use_grid_ylabels:
                    plot_results_kwargs_array[i][j]["ylabel"] = ""
                if i != 0 and use_grid_titles:
                    plot_results_kwargs_array[i][j]["title"] = ""
                if i != shape[0] - 1 and use_grid_xlabels:
                    plot_results_kwargs_array[i][j]["xlabel"] = ""
                if j != shape[1] - 1 and use_grid_hist_ylabels:
                    plot_results_kwargs_array[i][j]["hist_ylabel"] = ""

                # Plot results
                plot_results_kwargs = dict(
                    plot_results_kwargs_base, **plot_results_kwargs_array[i][j]
                )
                outpath_i = outpath.split(".")
                ext = outpath_i.pop(-1)
                outpath_i = ".".join(outpath_i)
                outpath_i = "".join([outpath_i, f"___arrloc_{i}_{j}.", ext])
                plot_results_kwargs["outpath"] = outpath_i
                plot_results(ax[i, j], **graph_array[i][j], **plot_results_kwargs)

    # Save figure
    if (
        "lumi_ratio" in plot_results_kwargs_base
        and plot_results_kwargs_base["lumi_ratio"] > 0.0
    ):
        outpath = outpath.replace(".pdf", "_rescaled.pdf")
    f.savefig(outpath)
