import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

try:
    from bokeh.plotting import ColumnDataSource, figure
    from bokeh.models import HoverTool, Range1d, FixedTicker
    from .bokeh_tools import FixedTickFormatter
    import pandas
    BOKEH_PANDAS_AVAILABLE = True
except ImportError:
    BOKEH_PANDAS_AVAILABLE = False


def updated_dict(dic1, dic2):
    if dic2 is not None:
        dic1.update(dic2)
    return dic1


class Band:

    def __init__(self, dna_size, migration_distance=None, ladder=None,
                 band_color="#000000", band_thickness=2, band_width=0.7,
                 label=None, label_fontdict=None, html=None):
        self.dna_size = dna_size
        if ladder is not None:
            migration_distance = ladder.dna_size_to_migration(dna_size)
        self.migration_distance = migration_distance
        self.band_color = band_color
        self.band_width = band_width
        self.band_thickness = band_thickness
        if label == "size":
            self.label = self.format_label(dna_size)
        else:
            self.label = label
        self.label_fontdict = label_fontdict
        self.html = html

    def format_label(self, dna_size):
        """Prettify the fragment size label, e.g. '13278' becomes '13.3k'"""
        if dna_size >= 1000:
            kilobases = np.round(dna_size / 1000.0, 1)
            number_fmt = "%d" if (kilobases == int(kilobases)) else "%.1f"
            return (number_fmt % kilobases) + "k"
        else:
            return "%d" % dna_size

    def plot_band(self, ax, x_coord):
        ax.plot([x_coord - self.band_width / 2.0,
                 x_coord + self.band_width / 2.0],
                [-self.migration_distance, -self.migration_distance],
                lw=self.band_thickness, c=self.band_color)

    def plot_label(self, ax, x_coord):
        if self.label is None:
            return
        fontdict = updated_dict({"color": 'white', 'family': 'sans-serif',
                                 "weight": "bold", "size": 10},
                                self.label_fontdict)
        ax.text(
            x_coord, -self.migration_distance, self.label,
            horizontalalignment="center",
            verticalalignment="center",
            fontdict=fontdict,
            transform=ax.transData,
            bbox=dict(boxstyle="round", fc=self.band_color, lw=0)
        )

    def plot(self, ax, x_coord):
        self.plot_band(ax, x_coord)
        self.plot_label(ax, x_coord)

    def to_json(self, ladder=None):
        return {prop: self.__dict__[prop] for prop in
                ["dna_size", "band_color", "band_width", "band_thickness",
                 "label", "html", "migration_distance"]}

    def modified(self, **attributes):
        new_obj = deepcopy(self)
        new_obj.__dict__.update(attributes)
        return new_obj


class BandsPattern:

    def __init__(self, bands, ladder=None, label=None, label_fontdict=None,
                 background_color=None, width=1.0, global_bands_props=None):
        self.bands = [
            Band(band, ladder=ladder) if isinstance(band, (int, float))
            else band
            for band in bands
        ]
        self.global_bands_props = ({} if global_bands_props is None
                                   else global_bands_props)
        self.label = label
        self.label_fontdict = label_fontdict
        self.background_color = background_color
        self.width = width
        self.initialize()

    def initialize(self):
        if len(self.bands) < 2:
            return
        self.dna_sizes, self.migration_distances = [
            np.array(e)
            for e in zip(*[
                (band.dna_size, band.migration_distance)
                for band in sorted(self.bands, key=lambda b: b.dna_size)
            ])
        ]
        self.migration_distance_span = (min(self.migration_distances),
                                        max(self.migration_distances))
        self.dna_size_span = min(self.dna_sizes), max(self.dna_sizes)

    def dna_size_to_migration(self, dna_sizes):
        return np.interp(dna_sizes, self.dna_sizes, self.migration_distances,
                         left=min(self.migration_distances),
                         right=max(self.migration_distances))

    def migration_to_dna_size(self, migration_distances):
        return np.interp(migration_distances, self.migration_distances[::-1],
                         self.dna_sizes[::-1], left=min(self.dna_sizes),
                         right=max(self.dna_sizes))

    def processed_bands(self):
        if self.global_bands_props == {}:
            return self.bands
        else:
            return [b.modified(**self.global_bands_props) for b in self.bands]

    def plot_background(self, ax, x_coord):
        if self.background_color is None:
            return
        ax.axvspan(xmin=x_coord - self.width / 2.0,
                   xmax=x_coord + self.width / 2.0,
                   color=self.background_color, zorder=-1000)
        # ax.fill_between([x_coord - self.width / 2.0,
        #                  x_coord + self.width / 2.0],
        #                 [-1000, -1000], facecolor=self.background_color,
        #                 linewidth=0, zorder=-1000)

    def plot_bands(self, ax, x_coord):
        for band in self.processed_bands():
            band.plot(ax, x_coord)

    def plot_label(self, ax, x_coord):
        if self.label in (None, ""):
            return
        fontdict = updated_dict({"color": "black", 'family': 'sans-serif',
                                 "weight": "bold", "size": 11, "rotation": 90},
                                self.label_fontdict)
        if 10 < fontdict["rotation"] < 80:
            alignment, shift = "left", -0.2
        else:
            alignment, shift = "center", 0
        ax.text(
            x_coord + shift, 0, "  " + self.label,
            horizontalalignment=alignment,
            verticalalignment="bottom", fontdict=fontdict,
            transform=ax.transData,
        )

    def plot(self, ax, x_coord):
        self.plot_background(ax, x_coord)
        self.plot_bands(ax, x_coord)
        self.plot_label(ax, x_coord)

    def merge_with(self, other):
        return self.modified(bands=self.bands + other.bands)

    def modified(self, **attributes):
        new_obj = deepcopy(self)
        new_obj.__dict__.update(attributes)
        new_obj.initialize()
        return new_obj

    def to_json(self):
        return


class BandsPatternsSet:

    def __init__(self, patterns, ladder=None, label=None, label_fontdict=None,
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

    def processed_patterns(self):
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

    def plot_ladder_ticks(self, ax):
        ticks = self.ladder_ticks
        if ticks is None:
            return
        if self.ladder is None:
            raise ValueError("Provide a `ladder` to BandsPatternsSet to enable"
                             " ladder ticks display.")
        fontdict = updated_dict(dict(size=7, rotation=90), self.ticks_fontdict)
        if isinstance(ticks, int):
            bmin, bmax = self.ladder.migration_distance_span
            migrations = np.linspace(bmin, bmax, ticks)
            ticks = [self.ladder.migration_to_dna_size(m) for m in migrations]
            ticks = [int(np.round(b, -2)) for b in ticks]  # round to 100
        yticks = [-self.ladder.dna_size_to_migration(b) for b in ticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(ticks, fontdict=fontdict)
        ax.yaxis.set_ticks_position('left')

    def plot_label(self, ax):
        if self.label is None:
            return
        fontdict = updated_dict({"weight": "bold", "size": 12},
                                self.label_fontdict)
        ax.set_ylabel(self.label, fontdict=fontdict, labelpad=8)

    def plot_patterns(self, ax):
        for i, pattern in enumerate(self.processed_patterns()):
            pattern.plot(ax, i + 1)

    def initialize_ax(self, ax):
        ax.set_frame_on(False)
        ax.set_yticks([])
        ax.set_xticks([])
        y1, y2 = self.ladder.migration_distance_span
        ax.set_ylim(-1.1 * y2, 0)

    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(0.5 * len(self.patterns), 3))
        self.initialize_ax(ax)
        self.plot_label(ax)
        self.plot_ladder_ticks(ax)
        self.plot_patterns(ax)
        return ax

    def plot_with_bokeh(self, visible_bands=12, band_width_pixels=40):

        if not BOKEH_PANDAS_AVAILABLE:
            raise ImportError("Install Bokeh and Pandas to use this feature")
        max_x = min(visible_bands, len(self.patterns) + 1)
        max_migration = self.ladder.migration_distances.max()
        mmin, mmax = self.ladder.migration_distance_span
        hw = 0.002 * abs(mmax - mmin)
        fig = figure(tools=[HoverTool(tooltips="@html", names=["bands"])] +
                           ["xwheel_zoom,xpan,reset, resize"],
                     plot_height=300, plot_width=band_width_pixels * max_x,
                     x_range=Range1d(0.5, max_x),  # labels,
                     y_range=Range1d(-1.1 * max_migration, 0),
                     logo=None, toolbar_location="right",
                     x_axis_location="above",
                     title_location="below",
                     title=self.label)
        fig.xaxis[0].ticker = FixedTicker(
            ticks=range(1, len(self.patterns) + 1))
        fig.xaxis[0].formatter = FixedTickFormatter(labels={
            i + 1: "" if (pattern.label is None) else pattern.label
            for i, pattern in enumerate(self.patterns)
        })
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
                for x_coord, pattern in enumerate(self.processed_patterns())
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
                for x_coord, pattern in enumerate(self.processed_patterns())
                for band in pattern.bands
            ])))

        fig.yaxis.visible = False
        fig.outline_line_color = None
        fig.grid.grid_line_color = None
        fig.xaxis.major_label_orientation = 0.6
        fig.axis.major_tick_in = 0
        fig.axis.major_tick_out = 2

        return fig
