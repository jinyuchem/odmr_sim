"""
Visualization utilities for ODMR simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Union, Dict


# Color schemes
COLORS = {
    'google': ['#4285F4', '#DB4437', '#F4B400', '#0F9D58', '#AB47BC',
               '#FF6F61', '#46BDC6', '#FF9E80', '#7E57C2', '#9CCC65'],
    'seven_level': ['#4285F4', '#DB4437', '#F4B400',  # GS states
                    '#4285F4', '#DB4437', '#F4B400',  # ES states (dotted)
                    '#0F9D58'],  # Singlet
    'readout': ['#AB47BC', '#46BDC6', '#FF9E80'],  # ES populations from different init
    'odmr': '#000080',  # Navy for ODMR spectrum
}

LINESTYLES_SEVEN_LEVEL = ['-', '-', '-', ':', ':', ':', '--']

LABELS_SEVEN_LEVEL = [
    "GS|0>", "GS|->", "GS|+>",
    "ES|0>", "ES|->", "ES|+>",
    "SS"
]


def setup_matplotlib(font_size: int = 10):
    """Set up matplotlib with consistent styling."""
    plt.rcParams.update({'font.size': font_size})


def plot_populations(
    t: np.ndarray,
    populations: np.ndarray,
    labels: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    linestyles: Optional[List[str]] = None,
    ax: Optional[plt.Axes] = None,
    xlabel: str = "Time (ns)",
    ylabel: str = "Population",
    title: Optional[str] = None,
    use_log_x: bool = True,
    xlim: Optional[tuple] = None,
    ylim: tuple = (-0.05, 1.15),
    legend_loc: str = 'upper left',
    legend_ncols: int = 2,
    time_unit: str = 'ns',
    show_legend: bool = True,
) -> plt.Axes:
    """Plot population dynamics over time.

    Parameters
    ----------
    t : np.ndarray
        Time array in seconds.
    populations : np.ndarray
        Population array with shape (n_times, n_states).
    labels : list of str, optional
        Labels for each state. If None and n_states=7, uses default labels.
    colors : list of str, optional
        Colors for each state.
    linestyles : list of str, optional
        Line styles for each state.
    ax : plt.Axes, optional
        Matplotlib axes to plot on. If None, creates new figure.
    xlabel, ylabel, title : str
        Axis labels and title.
    use_log_x : bool
        If True, use logarithmic x-axis.
    xlim, ylim : tuple
        Axis limits.
    legend_loc : str
        Legend location.
    legend_ncols : int
        Number of legend columns.
    time_unit : str
        Time unit for display ('ns', 'us', 's').
    show_legend : bool
        Whether to show the legend.

    Returns
    -------
    ax : plt.Axes
        The matplotlib axes.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    n_states = populations.shape[1]

    # Default labels for 7-level model
    if labels is None:
        if n_states == 7:
            labels = LABELS_SEVEN_LEVEL
        else:
            labels = [f"State {i}" for i in range(n_states)]

    # Default colors
    if colors is None:
        if n_states == 7:
            colors = COLORS['seven_level']
        else:
            colors = COLORS['google'][:n_states]

    # Default linestyles
    if linestyles is None:
        if n_states == 7:
            linestyles = LINESTYLES_SEVEN_LEVEL
        else:
            linestyles = ['-'] * n_states

    # Convert time units
    if time_unit == 'ns':
        t_plot = t * 1e9
    elif time_unit == 'us':
        t_plot = t * 1e6
    else:
        t_plot = t

    # Plot each state
    for j in range(n_states):
        ax.plot(
            t_plot, populations[:, j],
            color=colors[j % len(colors)],
            linestyle=linestyles[j % len(linestyles)],
            linewidth=1.5,
            label=labels[j]
        )

    # Formatting
    if use_log_x:
        ax.set_xscale('log')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title, fontsize=12)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    # Reference lines
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axhline(y=1, color='k', linewidth=0.5)

    # Tick formatting
    ax.tick_params(axis='both', direction='in')
    ax.tick_params(which='minor', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    if show_legend:
        ax.legend(
            fontsize=12, loc=legend_loc, edgecolor='black',
            borderpad=0.3, labelspacing=0.2, handlelength=1., handleheight=0.5,
            handletextpad=0.4, borderaxespad=0.3, columnspacing=0.6,
            ncols=legend_ncols
        )

    return ax


def plot_odmr_spectrum(
    frequencies: np.ndarray,
    contrast: np.ndarray,
    ax: Optional[plt.Axes] = None,
    xlabel: str = "MW Frequency (GHz)",
    ylabel: str = r"$\Delta I_{\mathrm{PL}} / I_{\mathrm{PL}}$ (%)",
    title: Optional[str] = None,
    color: str = '#000080',
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
    as_percentage: bool = True,
    negative_contrast: bool = True,  # Show as dip (experimental convention)
) -> plt.Axes:
    """Plot ODMR spectrum (contrast vs frequency).

    Parameters
    ----------
    frequencies : np.ndarray
        Microwave frequencies in GHz.
    contrast : np.ndarray
        ODMR contrast (fractional or percentage).
    ax : plt.Axes, optional
        Matplotlib axes.
    xlabel, ylabel, title : str
        Axis labels and title.
    color : str
        Line color.
    xlim, ylim : tuple
        Axis limits.
    as_percentage : bool
        If True, multiply contrast by 100 for display.
    negative_contrast : bool
        If True, negate contrast to show as dip (experimental convention).

    Returns
    -------
    ax : plt.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    contrast_plot = contrast * 100 if as_percentage else contrast
    if negative_contrast:
        contrast_plot = -contrast_plot

    ax.plot(frequencies, contrast_plot, color=color, linewidth=1.5)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title, fontsize=12)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.tick_params(axis='both', direction='in')
    ax.tick_params(which='minor', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    return ax


def plot_readout_comparison(
    results: Dict[str, dict],
    ax: Optional[plt.Axes] = None,
    xlabel: str = "Time (ns)",
    ylabel: str = r"$P_{\mathrm{ES, Total}}(t)$",
    title: Optional[str] = None,
    xlim: Optional[tuple] = None,
    ylim: Optional[tuple] = None,
) -> plt.Axes:
    """Plot readout comparison for different initial states.

    Parameters
    ----------
    results : dict
        Dictionary from ReadoutSimulation.run_comparison().
        Keys are initial state names, values are result dicts.
    ax : plt.Axes, optional
        Matplotlib axes.
    xlabel, ylabel, title : str
        Axis labels and title.
    xlim, ylim : tuple
        Axis limits.

    Returns
    -------
    ax : plt.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))

    colors = COLORS['readout']
    linestyles = ['-', '--', ':']
    labels = [
        r'Init. GS$|0\rangle$',
        r'Init. GS$|-\rangle$',
        r'Init. GS$|+\rangle$',
    ]

    for i, (state_name, result) in enumerate(results.items()):
        if result['es_total'] is not None:
            ax.plot(
                result['t_ns'],
                result['es_total'],
                color=colors[i % len(colors)],
                linestyle=linestyles[i % len(linestyles)],
                linewidth=1.5,
                label=labels[i]
            )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title, fontsize=12)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.tick_params(axis='both', direction='in')
    ax.tick_params(which='minor', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    ax.legend(
        fontsize=12, loc='best', edgecolor='black',
        borderpad=0.3, labelspacing=0.3, handlelength=1, handleheight=0.5,
        handletextpad=0.5, borderaxespad=0.3
    )

    return ax
