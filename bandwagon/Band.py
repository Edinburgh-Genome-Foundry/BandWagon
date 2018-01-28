import numpy as np
from copy import deepcopy
from .tools import updated_dict

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
