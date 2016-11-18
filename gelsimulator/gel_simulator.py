import numpy as np
from .tools import compute_digestion_bands


class GelSimulator:
    """A Gel plot generator.

    A gel simulator is defined by a `ladder` (a `GelLadder` object) using which
    it can plot bands on a matplotlib ax (e.g. bands produced by a digestion).

    Note that this was written as a class so that every little method can be
    overwritten, leaving you complete freedom to change the way the bands
    and the texts are rendered, simply by sub-classing GelSimulator with
    overwritten methods.

    Parameters
    ----------

    ladder

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
        ax.set_frame_on(False)
        ax.set_yticks([])
        ax.set_xticks([])
        y1, y2 = self.ladder.migration_distances_span()
        ax.set_ylim(-1.1 * y2, 0.1 * y2)

    def format_band_size_label(self, band_size):
        """Prettify the fragment size label, e.g. '13278' becomes '13.3k'"""
        if band_size >= 1000:
            kilobases = np.round(band_size / 1000.0, 1)
            number_fmt = "%d" if (kilobases == int(kilobases)) else "%.1f"
            return (number_fmt % kilobases) + "k"
        else:
            return "%d" % band_size

    def plot_band_rectangle(self, ax, band_size, x_coord, color="k",
                            width=0.6, band_thickness=4):
        """Plot one band (without label) on the matplotlib ax.

        At the moment it just draws a thick horizontal line at the specified
        location.

        Parameters
        ----------

        ax
          A matplotlib ax object

        band_size
          The size in basepairs of the DNA forming the band

        x_coord
          The x coordinate at which the pattern will be drawn.

        color
          Color of the band (can be any matplotlib-compatible color name)

        width

        band_thickness
          Thickness of the band in pixels. Beware that this means the band will
          have the same thickness undepending on the size

        """
        distance = self.ladder.compute_migration_distance(band_size)
        # patch = FancyBboxPatch([x_coord-0.4, -distance-0.2)
        ax.plot([x_coord - width/2.0, x_coord + width/2.0],
                [-distance, -distance],
                lw=band_thickness, c=color)

    def plot_band_size_label(self, ax, band_size, x_coord,
                             bg_color="black", label_fontdict=None):
        """Plot the band's label on the matplotlib ax.

        The label is meant to go over the result of `plot_band_rectangle`.

        Parameters
        ----------

        ax

        band_size
          Size in basepairs of the DNA forming the band

        x_coord
          x coordinate at which to the band is.

        label_fontdict
          A Matplotlib fontdict, e.g. {"color": "red", "size": 10}, to
          overwrite the label text defaults (white bold font)

        bg_color
          Color of the label's background

        """
        size_label = self.format_band_size_label(band_size)
        distance = self.ladder.compute_migration_distance(band_size)
        fontdict = {"color": 'white', 'family': 'sans-serif',
                    "weight": "bold", "size": 10}
        if label_fontdict is not None:
            fontdict.update(label_fontdict)
        ax.text(
            x_coord, -distance, size_label,
            horizontalalignment="center",
            verticalalignment="center",
            fontdict=fontdict,
            transform=ax.transData,
            bbox=dict(boxstyle="round", fc=bg_color, lw=0)
        )

    def plot_band(self, ax, band_size, x_coord, color="k",
                  width=0.3, with_size_labels=True, band_thickness=4,
                  label_fontdict=None):
        """Plot a band (horizontal line and size label) at the given location
        """
        self.plot_band_rectangle(ax, band_size, x_coord, color=color,
                                 width=width, band_thickness=band_thickness)
        if with_size_labels:
            self.plot_band_size_label(ax, band_size, x_coord, bg_color=color,
                                      label_fontdict=label_fontdict)

    def plot_horizontal_line(self, ax, band_size, xmin=0, xmax=1,
                             color="k", linewidth=1, linestyle="-", **kw):
        """Plot a horizontal line accross the ax at the given band size.

        Useful for comparing a series of band patterns with some expected
        band lengths.
        """
        y = -self.ladder.compute_migration_distance(band_size)
        ax.axhline(y, xmin=xmin, xmax=xmax, color=color, linewidth=linewidth,
                   linestyle=linestyle, **kw)

    def plot_band_pattern_label(self, ax, label, x_coord,
                                label_fontdict=None):
        """Write the label at the top of a band pattern
        """
        fontdict = {"color": "black", 'family': 'sans-serif', "weight": "bold",
                    "size": 13}
        if label_fontdict is not None:
            fontdict.update(label_fontdict)

        ax.text(
            x_coord, 0.05 * self.ladder.migration_distances_span()[1], label,
            horizontalalignment="center",
            verticalalignment="bottom",
            fontdict=fontdict,
            transform=ax.transData,
        )

    def color_band_zone(self, ax, x_coord, width=0.3, facecolor=None,
                        edgecolor="black", linewidth=0.5):
        ymax = self.ladder.migration_distances_span()[1]
        ax.fill_between([x_coord-width, x_coord+width], [-ymax, -ymax],
                        facecolor=facecolor, edgecolor=edgecolor,
                        linewidth=linewidth)

    def plot_bands_pattern(self, ax, bands, label=None, x_coord=1,
                           color="black", with_size_labels=True, width=0.3,
                           band_thickness=4, label_tilt_angle=0,
                           label_fontdict=None):

        """Plot all bands in a same column and add a label for the pattern
        """
        for band in bands:
            self.plot_band(ax=ax, band_size=band, x_coord=x_coord,
                           color=color, with_size_labels=with_size_labels,
                           width=width, band_thickness=band_thickness)

        if label is not None:
            self.plot_band_pattern_label(ax, label=label, x_coord=x_coord,
                                         label_fontdict=label_fontdict)



    def plot_ladder(self, ax, label="ladder", x_coord=1, color="r",
                    label_tilt_angle=0, with_size_labels=True,
                    label_fontdict=None):
        """Plot this simulator's ladder on a Matplotlib ax.

        (usually the first thing you want to plot)
        """
        bands = sorted(self.ladder.bands.keys())
        self.plot_bands_pattern(ax, bands, label=label, x_coord=x_coord,
                                color=color, with_size_labels=with_size_labels,
                                label_tilt_angle=label_tilt_angle,
                                label_fontdict=label_fontdict)

    def plot_digestion_result(self, ax, sequence, enzymes, linear=True,
                              label=None, x_coord=1, color="k",
                              with_size_labels=True, label_tilt_angle=0):
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
                                color=color, with_size_labels=with_size_labels,
                                label_tilt_angle=label_tilt_angle)
