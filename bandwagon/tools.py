"""Basic useful functions for computing bands and comparing patterns."""

from Bio.Seq import Seq
from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from Bio import Restriction
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from snapgene_reader import snapgene_file_to_seqrecord
import numpy as np


def random_dna_sequence(length):
    return "".join(["ATGC"[np.random.randint(0, 4)] for i in range(length)])


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(
        Seq(sequence, alphabet=DNAAlphabet()),
        id=id,
        name=name,
        features=list(features),
    )


def find_cut_sites(sequence, enzymes, linear=True):
    """Return the list of indices of cut sites in the given sequence."""
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    if hasattr(sequence, "seq"):
        sequence = sequence.seq
    batch = Restriction.RestrictionBatch(enzymes)
    matches = batch.search(sequence, linear=linear)
    return sorted(sum(matches.values(), []))


def compute_digestion_bands(sequence, enzymes, linear="auto_or_true"):
    """Return the band sizes [75, 2101, ...] resulting from enzymatic digestion

    Parameters
    ----------

    sequence
      Sequence to be digested. Either a string ("ATGC...") or a BioPython `Seq`

    enzymes
      list of all enzymes placed at the same time in the digestion mix
      e.g. `["EcoRI", "BamHI"]`

    linear
      True if the DNA fragment is linearized, False if it is circular. If you
      leave it to the default "auto_or_true", the function will read the
      topology from the sequence if it is a record with
      ``record.annotations['topology']`` set to "circular" or "linear", else it
      will default to true
    """
    if linear == "auto_or_true":
        linear = True
        if hasattr(sequence, "annotations"):
            if sequence.annotations.get("topology", None) == "circular":
                linear = False
    cut_sites = find_cut_sites(sequence, enzymes, linear=linear)
    cut_sites = [0] + cut_sites + [len(sequence)]
    bands_sizes = [
        end - start for (start, end) in zip(cut_sites, cut_sites[1:])
    ]
    if not linear and len(bands_sizes) > 1:
        bands_sizes[0] += bands_sizes.pop()
    return sorted(bands_sizes)


def annotate_record(
    seqrecord,
    location="full",
    feature_type="misc_feature",
    margin=0,
    **qualifiers
):
    """Add a feature to a Biopython SeqRecord. (also returns that same record)
    """
    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if (len(location) == 3) else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type,
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
        ax,
        iwidth,
        iheight,
        loc=10,  # means "center"
        bbox_to_anchor=[left, bottom, width, height],
        bbox_transform=ax.transData,
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



def load_record(
    record_file,
    topology="auto",
    default_topology="linear",
    id="auto",
    upperize=True,
    max_name_length=20,
    file_format=None,
):
    """Read a record (from many different input formats).

    Parameters
    ----------

    record_file
      A genbank file, a fasta file, a snapgene file, or a filelike object
      (at which case the format, genbank or fasta, must be given with
      ``file_format``)

    topology
      Either circular or linear or auto. If auto, then will attempt to read
      record.annotations['topology], and default to ``default_topology``.

    default_topology
      Default topology to use when topology is set to "auto" and a record
      has no designated topology.

    id
      Will be used for the record ID and name. If auto, the record id will
      be unchanged unless it is ".", " ", etc. at which case it will be
      replaced by the file name.

    upperize
      If true, the sequence will get upperized.

    max_name_length
      The name of the record will be truncated if too long to avoid Biopython
      exceptions being raised.

    file_format
      Indicates the file format for the parser, when record_file is a filelike
      object.
  
    """
    if file_format is not None:
        record = SeqIO.read(record_file, file_format)
    elif record_file.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(record_file, "genbank")
    elif record_file.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(record_file, "fasta")
    elif record_file.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(record_file)
    else:
        raise ValueError("Unknown format for file: %s" % record_file)
    if upperize:
        record = record.upper()
    if topology == "auto":
        set_record_topology(record, default_topology, pass_if_already_set=True)
    else:
        set_record_topology(record, topology)
    if id == "auto":
        id = record.id
        if id in [None, "", "<unknown id>", ".", " "]:
            id = os.path.splitext(os.path.basename(record_file))[0]
            record.name = id.replace(" ", "_")[:max_name_length]
        record.id = id
    elif id is not None:
        record.id = id
        record.name = id.replace(" ", "_")[:max_name_length]

    return record


def record_is_linear(record, default=True):
    """Return True if record.annotations['topology'] == 'linear'"""
    if "topology" not in record.annotations:
        return default
    else:
        return record.annotations["topology"] == "linear"

def set_record_topology(record, topology, pass_if_already_set=False):
    """Set record.annotations['topology'] (if not already set?)"""
    record_topology = record.annotations.get("topology", None)
    do_nothing = pass_if_already_set and (record_topology is not None)
    if not do_nothing:
        record.annotations["topology"] = topology