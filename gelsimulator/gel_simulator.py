from Bio import Restriction
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


def compute_digestion_bands(sequence, enzymes, linear=True):
    """Return the band sizes [75, 2101, ...] resulting from enzymatic digestion

    Parameters
    ----------

    sequence
      Sequence to be digested. Either a string ("ATGC...") or a BioPython `Seq`

    enzymes
      list of all enzymes placed at the same time in the digestion mix
      e.g. `["EcoRI", "BamHI"]`

    linear
      True if the DNA fragment is linearized, False if it is circular
    """
    if not isinstance(sequence, Seq):
        sequence = Seq(sequence)
    batch = Restriction.RestrictionBatch(enzymes)
    cut_sites = batch.search(sequence, linear=linear)
    cut_sites = [0] + sorted(sum(cut_sites.values(), [])) + [len(sequence)]
    bands_sizes = [end - start for (start, end)
                   in zip(cut_sites, cut_sites[1:])]
    if not linear:
        bands_sizes[0] += bands_sizes.pop()
    return sorted(bands_sizes)


class GelSimulator:
    """A Gel plot generator.

    A gel simulator is defined by a `ladder` (a `GelLadder` object) using which
    it can plot bands on a matplotlib ax (e.g. bands produced by a digestion).

    Note that this was written as a class so that every little method can be
    overwritten, leaving you complete freedom to change the way the bands
    and the texts are rendered, simply by sub-classing GelSimulator with
    overwritten methods.

    Examples
    --------

    >>> from gelsimulator import GelSimulator, LADDER_100_to_4k
    >>> import matplotlib.pyplot as plt
    >>> with open("example_sequence.txt", "r") as f:
    >>>     sequence = f.read()
    >>> enzymes = ["BamHI", "EcoRI", "EcoRV", "PstI", "SpeI", "XbaI"]
    >>> gel_simulator = GelSimulator(LADDER_100_to_4k)
    >>> fig, ax = plt.subplots(1, figsize=(1.1*(len(enzymes)+1), 5))
    >>> gel_simulator.format_ax(ax)
    >>> gel_simulator.plot_ladder(ax, x_coord=1)
    >>> for i, enzyme in enumerate(enzymes):
    >>>     gel_simulator.plot_digestion_result(ax, sequence, [enzyme],
    >>>                                         x_coord=i+2, label=enzyme)
    >>> fig.savefig("simple_gel_simulation.png", bbox_inches="tight")
    """

    def __init__(self, ladder):
        self.ladder = ladder

    def format_ax(self, ax):
        """Preformat the Matplotlib ax: remove the frame, set y-span."""
        ax.axis("off")
        y1, y2 = self.ladder.migration_distances_span
        ax.set_ylim(-1.1 * y2, 0.1 * y2)

    def format_fragment_size_label(self, band_size):
        """Prettify the fragment size label, e.g. '13278' becomes '13.3k'"""
        if band_size >= 1000:
            kilobases = np.round(band_size / 1000.0, 1)
            number_fmt = "%d" if (kilobases == int(kilobases)) else "%.1f"
            return (number_fmt % kilobases) + "k"
        else:
            return "%d" % band_size

    def plot_band_rectangle(self, ax, fragment_size, x_coord, color="k"):
        """Plot one band (without label) on the matplotlib ax.

        At the moment it just draws a thick horizontal line at the specified
        location
        """
        distance = self.ladder.compute_migration_distance(fragment_size)
        # patch = FancyBboxPatch([x_coord-0.4, -distance-0.2)
        ax.plot([x_coord - 0.3, x_coord + 0.3], [-distance, -distance],
                lw=4, c=color)

    def plot_band_size_label(self, ax, fragment_size, x_coord, color="k"):
        """Plot the band's label on the matplotlib ax.

        The label is meant to go over the result of `plot_band_rectangle`
        """
        size_label = self.format_fragment_size_label(fragment_size)
        distance = self.ladder.compute_migration_distance(fragment_size)
        ax.text(
            x_coord, -distance, size_label,
            horizontalalignment="center",
            verticalalignment="center",
            fontdict={"color": "white",
                      'family': 'sans-serif', "weight": "bold"},
            fontsize=10,
            transform=ax.transData,
            bbox=dict(boxstyle="round", fc=color, lw=0)
        )

    def plot_band(self, ax, fragment_size, x_coord, color="k",
                  with_size_label=True):
        """Plot a band (horizontal line and size label) at the given location
        """
        self.plot_band_rectangle(ax, fragment_size, x_coord, color=color)
        if with_size_label:
            self.plot_band_size_label(ax, fragment_size, x_coord, color=color)

    def plot_band_pattern_label(self, ax, label, x_coord):
        """Write the label at the top of a band pattern
        """
        ax.text(
            x_coord, 0.05 * self.ladder.migration_distances_span[1], label,
            horizontalalignment="center",
            verticalalignment="bottom",
            fontdict={"color": "black",
                      'family': 'sans-serif', "weight": "bold"},
            fontsize=13,
            transform=ax.transData,
        )

    def plot_bands_pattern(self, ax, bands, label=None, x_coord=1,
                           color="k"):

        """Plot all bands in a same column and add a label for the pattern
        """
        for band in bands:
            self.plot_band(ax=ax, fragment_size=band, x_coord=x_coord,
                           color=color)
        if label is not None:
            self.plot_band_pattern_label(ax, label=label, x_coord=x_coord)



    def plot_ladder(self, ax, label="ladder", x_coord=1, color="r"):
        """Plot this simulator's ladder on a Matplotlib ax.

        (usually the first thing you want to plot)
        """
        bands = sorted(self.ladder.bands.keys())
        self.plot_bands_pattern(ax, bands, label=label, x_coord=x_coord,
                                color=color)

    def plot_digestion_result(self, ax, sequence, enzymes, linear=True,
                              label=None, x_coord=1, color="k"):
        """ Plot the resulting band pattern of a digestion on a Matplotlib ax

        Parameters
        ----------

        ax
          A matplotlib ax

        sequence
          Either a string "ATGGTGC..." or a BioPython `Seq` object

        enzymes
          A list of restriction enzymes names  in the digestion mix,
          for instance `["EcoRI"]` or `["EcoRI", "BamHI"]`

        linear
          True if the DNA fragment is linear, False if it is circular

        label
          Label to be written at the top of the bands pattern

        x_coord
          x-coordinate of the column where the bands pattern will be drawn

        color
          Color of the bands and the labels.


        """
        bands = compute_digestion_bands(sequence, enzymes, linear=linear)
        self.plot_bands_pattern(ax, bands, label=label, x_coord=x_coord,
                                color=color)
