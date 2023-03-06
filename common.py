import matplotlib
import numpy as np
from cycler import cycler


def adjust_lightness(color, amount=0.5):
    """
    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib

    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def set_figure(fontsize=18, width=550.0, heightratio=None, height=None):
    r"""
    Parameters
    ----------
    fontsize : float
        sets the intended fontsize

    width : float
        sets the intended width in pts

    Notes
    -----
    To set equal to the columnwidth of the article:

    In the tex file 'Column width: \the\columnwidth' will print this size
    alternatively, '\message{Column width: \the\columnwidth}' will print to the log

    \linewidth should be used in place of \columnwidth if the figure is used
    within special enviroments (e.g. minipage)

    https://matplotlib.org/stable/tutorials/introductory/customizing.html
    https://scipy-cookbook.readthedocs.io/items/Matplotlib_LaTeX_Examples.html
    https://tex.stackexchange.com/questions/16942/difference-between-textwidth-linewidth-and-hsize
    """
    fig_width_pt = width
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches

    if heightratio is None:
        heightratio = (np.sqrt(5)-1.0)/2.0  # Aesthetic ratio
    if height is None:
        fig_height = fig_width*heightratio      # height in inches
    else:
        fig_height = height*inches_per_pt
    fig_size = [fig_width, fig_height]
    params = {'backend': 'pdf',
              'text.usetex': True,
              'text.latex.preamble': r"""
                                      \usepackage{amsmath}
                                      """,
              # fonts
              'font.family': 'serif',
              'font.serif': 'Nimbus Roman No9 L',
              # font sizes
              'axes.labelsize': fontsize,
              'font.size': fontsize,
              'axes.titlesize': fontsize,
              'legend.fontsize': fontsize,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              # figure size
              'figure.figsize': fig_size,
              'figure.constrained_layout.use': True,
              # line styling
              'lines.linewidth': 2,
              # legend
              'legend.frameon': False,
              # spines
              'axes.spines.top': True,
              'axes.spines.right': True,
              # saving
              'savefig.bbox': 'tight',
              'savefig.pad_inches': 1/72,
              # grid
              #'grid.color': '0.7',
              #'grid.linewidth': 0.1,
              # ticks
              'xtick.direction': 'inout',
              'ytick.direction': 'inout',
              'xtick.major.width':   0.6,
              'xtick.minor.width':   0.4,
              'ytick.major.width':   0.6,
              'ytick.minor.width':   0.4,
              'xtick.major.size':    3.5,
              'xtick.minor.size':    2,
              'ytick.major.size':    3.5,
              'ytick.minor.size':    2,
              # markers
              #'legend.numpoints': 2,
              # colors
              'axes.prop_cycle': cycler('color',
                                        ['tab:blue',
                                         'tab:red',
                                         'tab:green',
                                         'tab:orange',
                                         'tab:purple',
                                         'tab:brown',
                                         'tab:pink',
                                         'tab:gray',
                                         'tab:olive',
                                         'tab:cyan']),
              # axes
              'axes.linewidth': 0.6,
              }
    matplotlib.rcParams.update(params)
    from IPython.core.display import display, HTML
    display(HTML("<style>.container { width:80% !important; }</style>"))
    from IPython.display import set_matplotlib_formats
    set_matplotlib_formats('retina')
set_figure(width=500)
