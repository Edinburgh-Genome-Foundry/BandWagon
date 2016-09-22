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

    def __init__(self, ladder):
        self.ladder = ladder

    def format_ax(self, ax):
        ax.axis("off")
        y1, y2 = self.ladder.migration_distances_span
        ax.set_ylim(-1.1 * y2, 0.1 * y2)

    def format_fragment_size_label(self, band_size):
        if band_size >= 1000:
            kilobases = np.round(band_size / 1000.0, 1)
            number_fmt = "%d" if (kilobases == int(kilobases)) else "%.1f"
            return (number_fmt % kilobases) + "k"
        else:
            return "%d" % band_size

    def plot_band_rectangle(self, ax, fragment_size, x_coord, color="k"):
        """Bla bla bla"""

        distance = self.ladder.compute_migration_distance(fragment_size)
        # patch = FancyBboxPatch([x_coord-0.4, -distance-0.2)
        ax.plot([x_coord - 0.3, x_coord + 0.3], [-distance, -distance],
                lw=4, c=color)

    def plot_band_size_label(self, ax, fragment_size, x_coord, color="k"):
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
        self.plot_band_rectangle(ax, fragment_size, x_coord, color=color)
        if with_size_label:
            self.plot_band_size_label(ax, fragment_size, x_coord, color=color)

    def plot_bands_pattern(self, ax, bands, label=None, x_coord=1,
                           color="k"):
        for band in bands:
            self.plot_band(ax=ax, fragment_size=band, x_coord=x_coord,
                           color=color)
        if label is not None:
            self.plot_band_pattern_label(ax, label=label, x_coord=x_coord)

    def plot_band_pattern_label(self, ax, label, x_coord):
        ax.text(
            x_coord, 0.05 * self.ladder.migration_distances_span[1], label,
            horizontalalignment="center",
            verticalalignment="bottom",
            fontdict={"color": "black",
                      'family': 'sans-serif', "weight": "bold"},
            fontsize=13,
            transform=ax.transData,
        )

    def plot_ladder(self, ax, label="ladder", x_coord=1, color="r"):
        bands = sorted(self.ladder.bands.keys())
        self.plot_bands_pattern(ax, bands, label=label, x_coord=x_coord,
                                color=color)

    def plot_digestion_result(self, ax, sequence, enzymes, linear=True,
                              label=None, x_coord=1, color="k"):
        bands = compute_digestion_bands(sequence, enzymes, linear=linear)
        self.plot_bands_pattern(ax, bands, label=label, x_coord=x_coord,
                                color=color)
