"""Define Bandwagon's classes Band, BandPattern, BandsPatternsSet."""

import numpy as np
import matplotlib.pyplot as plt
from .tools import (updated_dict)

try:
    from bokeh.plotting import ColumnDataSource, figure
    from bokeh.models import HoverTool, Range1d, FixedTicker
    from bokeh.models import FuncTickFormatter
    import pandas
    BOKEH_PANDAS_AVAILABLE = True
except ImportError:
    BOKEH_PANDAS_AVAILABLE = False

from .BandsPattern import BandsPattern

class BandsPatternsSet:
    """ A set of band patterns, that will be plotted next to one another.

    Parameters
    ----------

    patterns
      A list of patterns, in the order in which they will be displayed. Each
      pattern can be either a BandPattern object or a list of Band objects or
      a list of DNA sizes.

    ladder
      A BandsPattern to use as a ladder.

    label
      A label that will be displayed vertically on the left of the final plot

    label_fontdict
      Font dictionnary for the label e.g. {'color': 'red', 'size': 7}

    global_patterns_props
      Properties that override that of each pattern in the set

    ladder_ticks
      List of dna sizes that will be represented by ticks on the left of the
      plot to guide the reading.

    ticks_fontdict
      Font dictionnary for the ticks indicating dna sizes.

    alternate_background_colors
      Colors that will be alternated for the backgrounds of the different
      columns. Leave to None for no background. If some patterns have a
      background color set this color will override the
      ``alternate_background_colors``.
    """

    def __init__(self, patterns, ladder, label=None, label_fontdict=None,
                 global_patterns_props=None, ladder_ticks=None,
                 ticks_fontdict=None,
                 alternate_background_colors=("#e2edff", "#fffae2")):
        self.patterns = [
            BandsPattern(p, ladder=ladder) if isinstance(p, (tuple, list))
            else p
            for p in patterns
        ]
        self.label = label
        self.label_fontdict = label_fontdict
        self.global_patterns_props = ({} if global_patterns_props is None
                                      else global_patterns_props)
        self.alternate_background_colors = alternate_background_colors
        self.ladder = ladder
        self.ladder_ticks = ladder_ticks
        self.ticks_fontdict = ticks_fontdict

    def _processed_patterns(self):
        """Versions of this set's patterns with attributes modified.

        The attributes are modified by ``global_patterns_props`` and by
        ``alternate_background_colors``.
        """
        new_patterns = []
        for i, pattern in enumerate(self.patterns):
            pattern = pattern.modified(**self.global_patterns_props)
            if pattern.background_color is None:
                if self.alternate_background_colors is not None:
                    ind = i % len(self.alternate_background_colors)
                    color = self.alternate_background_colors[ind]
                    pattern.background_color = color
            new_patterns.append(pattern)
        return new_patterns

    def _plot_ladder_ticks(self, ax):
        """Plot the ticks indicating the DNA sizes on the left of the plot"""
        ticks = self.ladder_ticks
        if ticks is None:
            return
        if self.ladder is None:
            raise ValueError("Provide a `ladder` to BandsPatternsSet to enable"
                             " ladder ticks display.")
        fontdict = updated_dict(
            dict(size=7, rotation=90, verticalalignment='center'),
            self.ticks_fontdict
        )
        if isinstance(ticks, int):
            bmin, bmax = self.ladder.migration_distances_span
            migrations = np.linspace(bmin, bmax, ticks)
            ticks = [self.ladder.migration_to_dna_size(m) for m in migrations]
            ticks = [int(np.round(b, -2)) for b in ticks]  # round to 100
        yticks = [-self.ladder.dna_size_to_migration(b) for b in ticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(ticks, fontdict=fontdict)
        ax.yaxis.set_ticks_position('left')

    def _plot_label(self, ax):
        """Plot the label of the figure on the left."""
        if self.label is None:
            return
        fontdict = updated_dict({"weight": "bold", "size": 12},
                                self.label_fontdict)
        ax.set_ylabel(self.label, fontdict=fontdict, labelpad=8)

    def _plot_patterns(self, ax):
        """Plot all band patterns side by side."""
        xmin, xmax = ax.get_xlim()
        patterns = self._processed_patterns()
        if xmax <= len(patterns) + 0.5:
            ax.set_xlim(xmax=len(patterns) + 0.5)
        for i, pattern in enumerate(self._processed_patterns()):
            pattern.plot(ax, i + 1)

    def _initialize_ax(self, ax):
        """Initialize the Matplotlib ax before plotting."""
        ax.set_frame_on(False)
        ax.set_yticks([])
        ax.set_xticks([])
        y1, y2 = self.ladder.migration_distances_span
        ax.set_ylim(-1.1 * y2, 0)

    def plot(self, ax=None):
        """Plot the band patterns on the given Matplotlib ax."""
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(0.5 * len(self.patterns), 3))
        self._initialize_ax(ax)
        self._plot_label(ax)
        self._plot_ladder_ticks(ax)
        self._plot_patterns(ax)
        ax.set_xlim(xmin=0.5)
        return ax

    def plot_with_bokeh(self, max_visible_patterns=12, band_width_pixels=40):
        """Return an interactive (browser-based) Bokeh figure of the patterns.

        Parameters
        ----------

        max_visible_patterns
          Max number of patterns that will be visible at the same time.
          A horizontal scroll will allow to see more patterns if there are more

        band_width_pixels
          Size of a band width in pixels on the screen. Said otherwise, the
          final figure will have a width of ``band_width_pixels * visible``
        """

        if not BOKEH_PANDAS_AVAILABLE:
            raise ImportError("Install Bokeh and Pandas to use this feature")
        max_x = min(max_visible_patterns, len(self.patterns) + 1)
        max_migration = self.ladder.migration_distances.max()
        mmin, mmax = self.ladder.migration_distances_span
        hw = 0.002 * abs(mmax - mmin)
        fig = figure(tools=[HoverTool(tooltips="@html", names=["bands"])] +
                           ["xwheel_zoom,xpan,reset"],
                     plot_height=300, plot_width=band_width_pixels * max_x,
                     x_range=Range1d(0.5, max_x),  # labels,
                     y_range=Range1d(-1.1 * max_migration, 0),
                     logo=None, toolbar_location="right",
                     x_axis_location="above",
                     title_location="below",
                     title=self.label)
        label_dict = {
            i + 1: "" if (pattern.label is None) else pattern.label
            for i, pattern in enumerate(self.patterns)
        }
        fig.xaxis[0].ticker = FixedTicker(
            ticks=list(range(1, len(self.patterns) + 1)))
        fig.xaxis.formatter = FuncTickFormatter(code="""
            var labels = %s;
            return labels[tick];
        """ % label_dict)
        fig.quad(
            name="backgrounds", top="top", bottom="bottom", left="left",
            right="right", color="color",
            source=ColumnDataSource(pandas.DataFrame.from_records([
                {
                    "left": x_coord + 1 - 0.5 * pattern.width,
                    "right": x_coord + 1 + 0.5 * pattern.width,
                    "top": 0,
                    "bottom": - 2 * max_migration,
                    "color": pattern.background_color,
                }
                for x_coord, pattern in enumerate(self._processed_patterns())
            ])))

        fig.quad(
            name="bands", top="top", bottom="bottom", left="left",
            right="right", color="color",
            source=ColumnDataSource(pandas.DataFrame.from_records([
                {
                    "left": np.round(x_coord + 1 - 0.5 * band.band_width, 2),
                    "right": np.round(x_coord + 1 + 0.5 * band.band_width, 2),
                    "top": - np.round(band.migration_distance -
                                      hw * band.band_thickness, 2),
                    "bottom": - np.round(band.migration_distance +
                                         hw * band.band_thickness, 2),
                    "color": band.band_color,
                    "html": band.html if band.html else (
                        band.label if band.label else (
                            "%d bp" % band.dna_size))
                }
                for x_coord, pattern in enumerate(self._processed_patterns())
                for band in pattern.bands
            ])))

        fig.yaxis.visible = False
        fig.outline_line_color = None
        fig.grid.grid_line_color = None
        fig.xaxis.major_label_orientation = 0.6
        fig.axis.major_tick_in = 0
        fig.axis.major_tick_out = 2

        return fig
