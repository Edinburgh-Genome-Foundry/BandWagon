"""Define Bandwagon's classes Band, BandPattern, BandsPatternsSet."""

import numpy as np
from scipy.interpolate import CubicSpline
from copy import deepcopy
from .tools import (place_inset_ax_in_data_coordinates, updated_dict)

from .Band import Band


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
        ``self.migration_distances``, ``self.migration_distances_span``,
        ``self.dna_size_span``, and indirectly the functions
        ``dna_size_to_migration`` and ``migration_to_dna_size``.
        """
        if len(self.bands) < 2:
            return
        self.dna_sizes, self.migration_distances = [
            np.array(e)
            for e in zip(*[
                (band.dna_size, band.migration_distance)
                for band in sorted(set(self.bands), key=lambda b: b.dna_size)
            ])
        ]
        self.migration_distances_span = (min(self.migration_distances),
                                         max(self.migration_distances))
        self.dna_size_span = min(self.dna_sizes), max(self.dna_sizes)
        self._dna_size_to_migration_interpolator = None
        self._migration_to_dna_size_interpolator = None

    def dna_size_to_migration(self, dna_sizes):
        """Return the migration distances for the given dna sizes"""
        if self._dna_size_to_migration_interpolator is None:
            self._dna_size_to_migration_interpolator = CubicSpline(
                self.dna_sizes, self.migration_distances,
                bc_type='natural'
            )
        return self._dna_size_to_migration_interpolator(dna_sizes)

    def migration_to_dna_size(self, migration_distances):
        """Return the dna sizes corresponding to the given migrations."""
        if self._migration_to_dna_size_interpolator is None:
            self._migration_to_dna_size_interpolator = CubicSpline(
                self.migration_distances[::-1], self.dna_sizes[::-1],
                bc_type='natural'
            )
        return self._migration_to_dna_size_interpolator(migration_distances)

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
            x_coord - self.width / 2.0, 0,
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
