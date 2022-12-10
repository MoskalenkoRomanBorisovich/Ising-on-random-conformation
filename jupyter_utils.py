import matplotlib
from IPython import get_ipython
import warnings
IPYTHON = get_ipython()

matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'

# beta symbol for matplotlib
BETA = r'$\beta$'

class matplotlib_pgf:
    """
    used to wrap matplotlib code to save figures as pgf

    Example:

    with matplotlib_pgf():
        plt.plot([1, 2, 3], [4, 5, 6])
        plt.grid(True)
        plt.savefig('test_plot.pgf')
    """
    def __inti__(self):
        self.RCDEFAULT = None
        
    def __enter__(self):
        self.RCDEFAULT = matplotlib.rcParamsDefault.copy()
        matplotlib.use('pgf')
        matplotlib.rcParams.update({
            "pgf.texsystem": "pdflatex",
            'font.family': 'serif',
            'text.usetex': True,
            'pgf.rcfonts': False,
        })

    def __exit__(self, exc_type, exc_value, traceback):
        global IPYTHON
        if self.RCDEFAULT is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                matplotlib.rcParams.update(self.RCDEFAULT)
            IPYTHON.magic("matplotlib inline")
            return True
        
        else:
            return False


def save_fig(fig: matplotlib.figure.Figure, fname: str):
    """saves figure as pgf and png

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        figure to save
    fname : str
        filename to save as
    """
    fig.savefig(f'{fname}.pgf', format='pgf', bbox_inches='tight')
    fig.savefig(f'{fname}.png', format='png', bbox_inches='tight')


def set_size(width_pt=345, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)