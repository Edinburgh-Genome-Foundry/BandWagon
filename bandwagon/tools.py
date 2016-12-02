"""Basic useful functions for computing bands and comparing patterns."""

from Bio.Seq import Seq
from Bio import Restriction
from mpl_toolkits.axes_grid.inset_locator import inset_axes

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
    if not linear and len(bands_sizes) > 1:
        bands_sizes[0] += bands_sizes.pop()
    return sorted(bands_sizes)


def bands_are_similar(b1, b2, tolerance):
    return float(abs(b1 - b2)) / min(b1, b2) < tolerance


def merge_bands_in_pattern(bands, tolerance):
    bands = sorted(bands)
    new_bands = [bands[0]]
    for band in bands[1:]:
        if bands_are_similar(band, new_bands[-1], tolerance):
            new_bands[-1] = 0.5 * (band + new_bands[-1])
        else:
            new_bands.append(band)
    return new_bands


def bands_patterns_are_similar(bands1, bands2,
                               tolerance=0.3, merge_tolerance=0.2):
    if (bands1 == []) or (bands2 == []):
        return False
    bands1 = merge_bands_in_pattern(bands1, tolerance=merge_tolerance)
    bands2 = merge_bands_in_pattern(bands2, tolerance=merge_tolerance)
    return (len(bands1) == len(bands2)) and all(
        bands_are_similar(b1, b2, tolerance)
        for b1, b2 in zip(sorted(bands1), sorted(bands2))
    )

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
