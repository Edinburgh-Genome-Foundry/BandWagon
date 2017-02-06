"""Define Bandwagon's classes Band, BandPattern, BandsPatternsSet."""

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from .tools import place_inset_ax_in_data_coordinates, updated_dict

try:
    from bokeh.plotting import ColumnDataSource, figure
    from bokeh.models import HoverTool, Range1d, FixedTicker
    from .bokeh_tools import FixedTickFormatter
    import pandas
    BOKEH_PANDAS_AVAILABLE = True
except ImportError:
    BOKEH_PANDAS_AVAILABLE = False





class Band:
    """A band in a migration pattern.

    This class controls all visual aspects of the band (form, color, label...).

    Parameters
    ----------

    dna_size

    migration_distance
      Observed migration distance of the band. Useful if you are using this
      band to define a ladder, otherwise keep to None.

    ladder
      Ladder used to determine the migration distance of the band, if no
      ``migration_distance`` was provided.

    band_color
      Color of the band. Can be a color name (like 'black', 'red'...) or
      an hexadecimal code like '#001aff' (this format is mandatory if you
      want to export to browser-based display.

    band_thickness
      Thickness of the band in pixels.

    band_width
      Proportion of the column width that is occupied by the band (the columns
      for each band pattern have a width of 1.0)

    label
      Label to be printed on the band, if any provided. Set to '=size' for
      printing the bands sizes

    label_fontdict
      A dict indicating the format of the label if any provided. For isinstance
      ``{'size': 7, 'weight':'bold', 'color': '#0011aa'}``

    html
      Html appearing when hovering the band.
    """
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
        self.label = label
        self.label_fontdict = label_fontdict
        self.html = html
        self.initialize()

    def initialize(self):
        if self.label == "=size":
            self.label = self._format_dna_size(self.dna_size)

    @staticmethod
    def _format_dna_size(dna_size):
        """Prettify the fragment size label, e.g. '13278' becomes '13.3k'"""
        if dna_size >= 1000:
            kilobases = np.round(dna_size / 1000.0, 1)
            number_fmt = "%d" if (kilobases == int(kilobases)) else "%.1f"
            return (number_fmt % kilobases) + "k"
        else:
            return "%d" % dna_size

    def _plot_band(self, ax, x_coord):
        """Plot the band's line on the matplotlib ax at the x-coordinate."""
        ax.plot([x_coord - self.band_width / 2.0,
                 x_coord + self.band_width / 2.0],
                [-self.migration_distance, -self.migration_distance],
                lw=self.band_thickness, c=self.band_color)

    def _plot_label(self, ax, x_coord):
        """Plot the label on (top of) the band."""
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
        """Plot the band's line and label on the ax at the x-coordinate."""
        self._plot_band(ax, x_coord)
        self._plot_label(ax, x_coord)

    def to_json(self, ladder=None):
        """Return a dictionnary version of the band (for bakeh plotting)"""
        return {prop: self.__dict__[prop] for prop in
                ["dna_size", "band_color", "band_width", "band_thickness",
                 "label", "html", "migration_distance"]}

    def modified(self, **attributes):
        """Return a new version of this band, with modified attributes."""
        new_obj = deepcopy(self)
        new_obj.__dict__.update(attributes)
        new_obj.initialize()
        return new_obj


class BandsPattern:
    """Bands set forming a migration pattern. Also used to define ladders.

    Parameters
    ----------

    bands
      Either a list of ``Band`` objects, or a list of DNA sizes e.g.
      ``[800, 654, 1752]``. If the latter case, a ``ladder`` must be provided
      to determine the elements

    ladder
      A BandsPattern to use as a ladder in case ``bands`` is a list of DNA
      sizes. Else leave it to None.

    label
      Label to be printed on top of the band pattern

    label_fontdict
      A dict indicating the format of the label if any provided. For instance
      ``{'size': 7, 'rotation':30, 'color': '#0011aa'}``

    corner_note
      Note to be printed in small prints in the corner of the band lane.
      Useful to print e.g. the total sum of all bands, or some remark.

    corner_note_fontdict
      A dict indicating the format of the label if any provided. For instance
      ``{'size': 6, 'color': '#0011aa'}``

    background_color
      Background color of the column in which the pattern is plotted. Either
      the name of a color or a html code ('#0012e4').

    width
      Width of the column (better keep to 1.0 if you want my humble opinion)

    global_bands_props
      Dictionnary of band properties overwriting that of all bands in the
      pattern, e.g. ``{'color': '#012ff45'}``

    gel_image
      A numpy array of size w x h and values in 0-1 representing a gel image.
      This image will be superimposed on the right of the pattern's color.

    gel_image_width
      Width of the gel_image display (remember that the width of the column is
      1.0, so 0.2 means this 'photo' will occupy 1/5 of the column, the rest
      being occupied by the pattern plot)
     """
    def __init__(self, bands, ladder=None, label=None, label_fontdict=None,
                 corner_note=None,  corner_note_fontdict=None,
                 background_color=None, width=1.0, global_bands_props=None,
                 gel_image=None, gel_image_width=0.2):
        self.bands = [
            Band(band, ladder=ladder) if isinstance(band, (int, float))
            else band
            for band in bands
        ]
        self.global_bands_props = ({} if global_bands_props is None
                                   else global_bands_props)
        self.ladder = ladder
        self.label = label
        self.label_fontdict = label_fontdict
        self.background_color = background_color
        self.gel_image = gel_image
        self.gel_image_width = gel_image_width
        self.width = width
        self.corner_note = corner_note
        self.corner_note_fontdict = corner_note_fontdict
        self.initialize()

    def initialize(self):
        """Create the variables necessary to compute some properties.

        Concerned properties are ``self.dna_sizes``,
        ``self.migration_distances``, ``self.migration_distance_span``,
        ``self.dna_size_span``, and indirectly the functions
        ``dna_size_to_migration`` and ``migration_to_dna_size``.
        """
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
        """Return the migration distances for the given dna sizes"""
        return np.interp(dna_sizes, self.dna_sizes, self.migration_distances,
                         left=max(self.migration_distances),
                         right=min(self.migration_distances))

    def migration_to_dna_size(self, migration_distances):
        """Return the dna sizes corresponding to the given migrations."""
        return np.interp(migration_distances, self.migration_distances[::-1],
                         self.dna_sizes[::-1], left=min(self.dna_sizes),
                         right=max(self.dna_sizes))

    def _processed_bands(self):
        """Return the bands modified by ``self.global_bands_props``"""
        if self.global_bands_props == {}:
            return self.bands
        else:
            return [b.modified(**self.global_bands_props) for b in self.bands]

    def _plot_background(self, ax, x_coord):
        """Return the bands modified by ``self.global_bands_props``"""
        if self.background_color is None:
            return
        ax.axvspan(xmin=x_coord - self.width / 2.0,
                   xmax=x_coord + self.width / 2.0,
                   color=self.background_color, zorder=-1000)

    def _plot_gel_image(self, ax, x_coord):
        """Plot the gel image ('photo') at the given coordinates"""
        if self.gel_image is None:
            return
        # bottom left width height
        bbox = (x_coord + 0.45 - self.gel_image_width,
                0.95 * ax.get_ylim()[0],
                self.gel_image_width,
                0.90 * abs(ax.get_ylim()[0]))
        new_ax = place_inset_ax_in_data_coordinates(ax, bbox=bbox)
        # new_ax.plot([0,1,2,3], [0,1,0,1])
        new_ax.axis("off")
        new_ax.imshow(+self.gel_image, interpolation="none", aspect='auto')

    def _plot_bands(self, ax, x_coord):
        """Plot all bands at the given x_coordinate on the ax."""
        for band in self._processed_bands():
            band.plot(ax, x_coord)

    def _plot_label(self, ax, x_coord):
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

    def _plot_corner_note(self, ax, x_coord):
        if self.corner_note in (None, ""):
            return
        fontdict = updated_dict({"size": 5.5, "rotation": 90},
                                self.corner_note_fontdict)
        ax.text(
            x_coord - self.width/2.0, 0,
            self.corner_note + " ", horizontalalignment="left",
            verticalalignment="top", fontdict=fontdict,
            transform=ax.transData,
        )

    def plot(self, ax, x_coord):
        """Plot background, bands, label, gel_image at the given x_coord."""
        self._plot_background(ax, x_coord)
        self._plot_bands(ax, x_coord)
        self._plot_label(ax, x_coord)
        self._plot_gel_image(ax, x_coord)
        self._plot_corner_note(ax, x_coord)



    def merge_with(self, other):
        """Merge this band pattern with another pattern's bands"""
        return self.modified(bands=self.bands + other.bands)

    def modified(self, **attributes):
        """Return a version of this bands pattern with modified attributes."""
        new_obj = deepcopy(self)
        new_obj.__dict__.update(attributes)
        new_obj.initialize()
        return new_obj


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

    def _plot_label(self, ax):
        """Plot the label of the figure on the left."""
        if self.label is None:
            return
        fontdict = updated_dict({"weight": "bold", "size": 12},
                                self.label_fontdict)
        ax.set_ylabel(self.label, fontdict=fontdict, labelpad=8)

    def _plot_patterns(self, ax):
        """Plot all band patterns side by side."""
        for i, pattern in enumerate(self._processed_patterns()):
            pattern.plot(ax, i + 1)

    def _initialize_ax(self, ax):
        """Initialize the Matplotlib ax before plotting."""
        ax.set_frame_on(False)
        ax.set_yticks([])
        ax.set_xticks([])
        y1, y2 = self.ladder.migration_distance_span
        ax.set_ylim(-1.1 * y2, 0)

    def plot(self, ax=None):
        """Plot the band patterns on the given Matplotlib ax."""
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(0.5 * len(self.patterns), 3))
        self._initialize_ax(ax)
        self._plot_label(ax)
        self._plot_ladder_ticks(ax)
        self._plot_patterns(ax)
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
