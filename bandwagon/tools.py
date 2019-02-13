"""Basic useful functions for computing bands and comparing patterns."""

from Bio.Seq import Seq
from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import numpy as np

def random_dna_sequence(length):
    return "".join(["ATGC"[np.random.randint(0, 4)] for i in range(length)])

def sequence_to_biopython_record(sequence, id='<unknown id>',
                                 name='<unknown name>', features=()):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(Seq(sequence, alphabet=DNAAlphabet()),
                     id=id, name=name, features=list(features))

def find_cut_sites(sequence, enzymes, linear=True):
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    if hasattr(sequence, 'seq'):
        sequence = sequence.seq
    batch = Restriction.RestrictionBatch(enzymes)
    matches = batch.search(sequence, linear=linear)
    return sorted(sum(matches.values(), []))

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
    cut_sites = find_cut_sites(sequence, enzymes, linear=linear)
    cut_sites = [0] + cut_sites + [len(sequence)]
    bands_sizes = [end - start for (start, end)
                   in zip(cut_sites, cut_sites[1:])]
    if not linear and len(bands_sizes) > 1:
        bands_sizes[0] += bands_sizes.pop()
    return sorted(bands_sizes)

def annotate_record(seqrecord, location="full", feature_type="misc_feature",
                margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord. (also returns that same record)
    """
    if location == "full":
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if (len(location) == 3) else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )
    return seqrecord


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

def load_record(filename, linear=True, name="unnamed", fmt='auto'):
    """Load a FASTA/Genbank/... record"""
    if fmt is not 'auto':
        record = SeqIO.read(filename, fmt)
    elif filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(('fa', 'fasta')):
        record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError('Unknown format for file: %s' % filename)
    record.linear = linear
    if name != "unnamed":
        record.id = name
        record.name = name.replace(" ", "_")[:20]
    return record