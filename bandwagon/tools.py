"""Basic useful functions for computing bands and comparing patterns."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import numpy as np

def np_spline_interpolator(xx, yy, npoints=3, deg=2):
    """A poorman's spline interpolator, not using Scipy.

    It is slower and more naive than Scipy's actual interpolator, but really
    does its job here of interpolating ladders.
    """
    radius = int((npoints - 1) / 2)
    centers = xx[radius: -radius]
    print (centers)
    interpolators = [
        np.polyfit(xx[i-radius: i+radius+1], yy[i-radius: i+radius+1], 2)
        for i in range(radius, len(xx) - radius)
    ]
    def f(x):
        ii = np.interp(np.array(x), centers, range(len(centers)),
                       left=0, right=(len(centers)-1))
        if hasattr(ii, '__iter__'):
            return np.array([
                np.polyval(interpolators[i], _x)
                for (i, _x) in zip(ii.astype(int), x)
            ])
        else:
            return np.polyval(interpolators[int(ii)], x)
    return f

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
    if isinstance(sequence, SeqRecord):
        sequence = sequence.seq
    if not isinstance(sequence, Seq):
        sequence = Seq(sequence)

    batch = Restriction.RestrictionBatch(enzymes)
    cut_sites = batch.search(sequence, linear=linear)
    cut_sites = [0] + sorted(sum(cut_sites.values(), [])) + [len(sequence)]
    bands_sizes = [end - start for (start, end)
                   in zip(cut_sites, cut_sites[1:])]
    if not linear and len(bands_sizes) > 1:
        bands_sizes[0] += bands_sizes.pop()
    return sorted(bands_sizes)


def place_inset_ax_in_data_coordinates(ax, bbox):
    """Return an ax inset in the given ax at the given bbox in
    data coordinates (left, bottom, width, height)"""

    left, bottom, width, height = bbox
    pixels_data_00 = ax.transData.transform([0, 0])
    pixels_data_wh = ax.transData.transform([width, height])
    iwidth, iheight = (pixels_data_wh - pixels_data_00) / ax.figure.dpi
    return inset_axes(
        ax, iwidth, iheight,
        loc=10,  # means "center"
        bbox_to_anchor=[left, bottom, width, height],
        bbox_transform=ax.transData
    )

def updated_dict(dic1, dic2):
    """Return dic1 updated with dic2 if dic2 is not None.

    Example
    -------

    >>> my_dict = updated_dict({"size": 7, "color": "r"}, some_defaults_dict)
    """
    if dic2 is not None:
        dic1.update(dic2)
    return dic1
