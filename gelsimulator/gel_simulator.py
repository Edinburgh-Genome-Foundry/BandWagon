from collections import OrderedDict
import numpy as np
from .tools import compute_digestion_bands
import matplotlib.pyplot as plt


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
        distance = self.ladder.band_size_to_migration(band_size)
        ax.plot([x_coord - width / 2.0, x_coord + width / 2.0],
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
        distance = self.ladder.band_size_to_migration(band_size)
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
        y = -self.ladder.band_size_to_migration(band_size)
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
                        edgecolor="black", linewidth=0.5, zorder=-1000,
                        ladder_only=False):
        ymax = self.ladder.migration_distances_span()[1]
        if not ladder_only:
            ymax *= 3
        ax.fill_between([x_coord - width, x_coord + width], [-ymax, -ymax],
                        facecolor=facecolor, edgecolor=edgecolor,
                        linewidth=linewidth, zorder=zorder)

    def plot_bands_pattern(self, ax, bands, label=None, x_coord=1,
                           color="black", with_size_labels=True, width=0.3,
                           band_thickness=2, label_tilt_angle=0,
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

    def plot_bands_patterns(self, patterns, ax=None, color='black',
                            with_size_labels=False, band_thickness=2,
                            plot_ladder=False, ladder_label="L",
                            plot_ladder_ticks=False,
                            plot_labels=True, patterns_labels_fontdict=None,
                            ylabel=None, ylabel_fontdict=None,
                            background_colors=("#e2edff", "#fffae2")):
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(0.5 * len(patterns), 3))
        if not isinstance(patterns, list):
            sort = not isinstance(patterns, OrderedDict)
            patterns = list(patterns.keys())
            if sort:
                patterns = sorted(patterns)
        self.format_ax(ax)
        if ylabel is not None:
            yfontdict = dict(weight="bold", size=10)
            if ylabel_fontdict is not None:
                yfontdict.update(ylabel_fontdict)
            ax.set_ylabel(ylabel, fontdict=yfontdict, labelpad=8)
        if plot_ladder:
            self.plot_ladder(ax, label=ladder_label,
                             label_fontdict=dict(size=7, rotation=20),
                             with_size_labels=with_size_labels,
                             band_thickness=band_thickness)
        if plot_ladder_ticks:
            self.plot_ladder_ticks(ax)
        for i, (name, pattern) in enumerate(patterns):
            x_coord = i + (2 if plot_ladder else 1)
            bgcolor = background_colors[i % 2]
            self.color_band_zone(ax, x_coord=x_coord, facecolor=bgcolor,
                                 linewidth=0, width=.5, zorder=-10000)
            pfontdict = dict(size=7, rotation=20)
            if patterns_labels_fontdict is not None:
                pfontdict.update(patterns_labels_fontdict)
            self.plot_bands_pattern(ax, pattern, x_coord=x_coord,
                                    label=name if plot_labels else None,
                                    label_fontdict=pfontdict,
                                    with_size_labels=with_size_labels,
                                    band_thickness=band_thickness)
        return ax

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

    def plot_ladder_ticks(self, ax, bands=3, fontdict=None):
        """

        Parameters
        ----------

        bands
          can be a list or a number

        """
        ticks_fontdict = dict(size=6, rotation=90)
        if fontdict is not None:
            ticks_fontdict.update(fontdict)
        if isinstance(bands, int):
            bmin, bmax = self.ladder.migration_distances_span()
            migrations = np.linspace(bmin, bmax, bands)
            bands = [self.ladder.migration_to_band_size(m) for m in migrations]
            bands = [int(np.round(b, -2)) for b in bands]  # round to 100
        yticks = [-self.ladder.band_size_to_migration(b) for b in bands]
        ax.set_yticks(yticks)
        ax.set_yticklabels(bands, fontdict=ticks_fontdict)
        ax.yaxis.set_ticks_position('left')


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

        def plot_digestions_results(self, sequences, enzymes, linear=True,
                                    ylabel="auto",
                                    **plot_bands_patterns_kwargs):
            """ bla.

            Parameters
            ----------

            """
            if not isinstance(sequences, OrderedDict):
                sequences = OrderedDict([e for e in sorted(sequences.items())])
            patterns = OrderedDict([
                (name, compute_digestion_bands(sequence, enzymes,
                                               linear=linear))
                for name, sequence in sequences.items()
            ])
            if plot_bands_patterns_kwargs.get(ylabel, None) is None:
                plot_bands_patterns_kwargs["ylabel"] = " + ".join(enzymes)
            ax = self.plot_bands_patterns(patterns,
                                          **plot_bands_patterns_kwargs)
            return ax
