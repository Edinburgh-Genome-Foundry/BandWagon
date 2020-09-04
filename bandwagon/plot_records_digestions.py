from copy import deepcopy
from collections import OrderedDict
import itertools

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

from .Band import Band
from .BandsPattern import BandsPattern
from .BandsPatternsSet import BandsPatternsSet
from .tools import (
    find_cut_sites,
    annotate_record,
    compute_digestion_bands,
    record_is_linear,
)

try:
    from dna_features_viewer import BiopythonTranslator

    class DigestGraphicTranslator(BiopythonTranslator):
        def compute_feature_color(self, feature):
            if "band_label" in feature.qualifiers:
                return "#e5d740"
            elif "source" in feature.qualifiers:
                return "#ff8282"
            else:
                return "#eaedff"


except ImportError:

    class DigestGraphicTranslator:
        """Not available, install dna_features_viewer.

        pip install dna_features_viewer
        """

        def __init__(self, *args, **kwargs):
            raise ImportError(
                "You need to install dna_features_viewer to plot"
                " annotated digestion graphics. Use "
                "'pip install dna_features_viewer'"
            )


def annotate_digestion_bands(record, enzymes, ladder):
    """Annotate the record to indicate the regions corresponding to bands.
    """

    linear = record_is_linear(record, default=False)
    sequence = record.seq
    all_cuts = sorted(
        [0, len(record)] + find_cut_sites(sequence, enzymes, linear=linear)
    )
    bands = list(zip(all_cuts, all_cuts[1:]))
    if (not linear) and len(bands) > 1:
        start, end = bands.pop()
        band0 = [-(end - start), bands[0][1]]
        if bands == []:
            bands = [band0]
        else:
            bands[0] = band0
    sorted_bands = sorted(bands, key=lambda b: b[0] - b[1])
    new_record = deepcopy(record)
    for (band, label) in zip(sorted_bands, "abcdefghijkl"):
        band_size = abs(band[1] - band[0])
        formatted_size = Band._format_dna_size(band_size)
        annotate_record(
            new_record,
            location=band,
            label="%s - %s" % (label, formatted_size),
            feature_type="misc_feature",
            band_label=label,
            band_size=band_size,
        )
    return new_record


def plot_record_digestion(record_digestion, ladder, record_label, digestion_label):
    """Plot the digestion along with a schema of cuts locations in the record.

    Parameters
    ----------
    record_digestion
      Biopython record with features indicating bands.

    ladder

    record_label
      Label to use as the title of the record digestion plot.

    digestion_label
      Label to use as the title of the bands pattern(s).

    """
    gs = gridspec.GridSpec(4, 10)

    ax_bands, ax_record, ax_digest = [
        plt.subplot(z) for z in (gs[:, 0], gs[:3, 3:], gs[3, 3:])
    ]
    ax_digest.figure.set_size_inches(20, 7)
    ax_digest.set_title("BANDS")
    topology = "linear" if record_is_linear(record_digestion) else "circular"
    ax_record.set_title("%s (%s)" % (record_label, topology), fontsize=22)

    bands = sorted(
        [
            (
                feature.qualifiers["band_label"],
                feature.qualifiers["band_size"],
                feature.location.start,
            )
            for feature in record_digestion.features
            if feature.qualifiers.get("band_label", False)
        ]
    )

    for _, size, start in bands:
        for ax in (ax_digest, ax_record):
            ax.axvline(start, ls=":", color="k", lw=0.5)
            ax.axvline(start + size, ls=":", color="k", lw=0.5)

    gr_record, gr_digestion = [
        DigestGraphicTranslator([fl]).translate_record(record_digestion)
        for fl in [
            lambda f: ("band_label" not in f.qualifiers)
            and ("homology" not in f.qualifiers.get("label", [])),
            lambda f: ("band_label" in f.qualifiers),
        ]
    ]
    gr_digestion.split_overflowing_features_circularly()
    gr_digestion.plot(ax=ax_digest)
    gr_record.plot(ax=ax_record, with_ruler=False)
    linear = record_is_linear(record_digestion, default=False)
    pattern = BandsPattern(
        [Band(dnasize, ladder=ladder, label=label) for label, dnasize, _ in bands],
        ladder=ladder,
        topology="linear" if linear else "circular",
    )
    patternset = BandsPatternsSet(
        [pattern],
        ladder=ladder,
        ladder_ticks=3,
        ticks_fontdict=dict(size=12),
        label=digestion_label,
        label_fontdict=dict(size=18),
    )
    patternset.plot(ax_bands)
    return (ax_bands, ax_record, ax_digest)


def plot_records_digestions(
    target, ladder, records_and_digestions=None, records=None, digestions=None
):
    """Plot records digestions in a multipage PDF file.

    Parameters
    ----------

    records_and_digestions
      A list [(record, digestion), ...] where record is a Genbank record and
      digestions a tuple of enzymes. Alternatively, a list of ``records`` and
      ``digestions`` can be provided.

    records
      List of biopython records whose digests need to be digested by the
      ``digestions``.

    digestions
      List of tuples of enzyme names e.g. ``[('BsaI',),`('MfeI', 'BsmBI')...]``
      representing the digestions for the provided ``records``.

    ladder
      A Ladder object representing the ladder to be used for placing the bands.

    target
      path to a PDF file, or file-like object.

    full_report
    """
    if records_and_digestions is None:
        records_and_digestions = itertools.product(records, digestions)
    annotated_records = OrderedDict()
    with PdfPages(target) as details_pdf:
        for record, enzymes in records_and_digestions:
            record_label = record.id
            if record_label not in annotated_records:
                annotated_records[record_label] = OrderedDict()
            enzymes_label = " + ".join(sorted(enzymes))
            basename = "%s--%s" % (record_label.replace(" ", "_"), "+".join(enzymes),)
            record_digestion = annotate_digestion_bands(record, enzymes, ladder)
            record_digestion.id = basename
            annotated_records[record_label][enzymes_label] = record_digestion
            (ax, _, _) = plot_record_digestion(
                record_digestion, ladder, record_label, enzymes_label
            )
            details_pdf.savefig(ax.figure, bbox_inches="tight")
            plt.close(ax.figure)
    return annotated_records


def plot_all_digestion_patterns(
    records,
    digestions,
    ladder,
    axes=None,
    group_by="digestions",
    show_band_sizes=False,
    plot_ladder=False,
):
    """Plot a grid (RECORD x DIGESTION) of predicted records digestions.

    Parameters
    ----------
    records
      List of biopython records whose digests need to be digested by the
      ``digestions``.

    digestions
      List of tuples of enzyme names e.g. ``[('BsaI',),`('MfeI', 'BsmBI')...]``
      representing the digestions for the provided ``records``.

    ladder
      A Ladder object representing the ladder to be used for placing the bands
    axes.
      If None, new axes will be created.

    group_by
      Either "sequence" or "digestion". Determines the plotting order.

    show_band_sizes
      If true, the band sizes will be printed on each band.

    """
    all_patterns = OrderedDict()
    for record in records:
        topology = record.annotations.get("topology", "linear")
        linear = topology == "linear"
        for enzymes in digestions:
            enzymes_label = " + ".join(sorted(enzymes))
            bands = compute_digestion_bands(record.seq, enzymes, linear=linear)
            cuts = find_cut_sites(record.seq, enzymes)
            bands = sorted(bands)
            if group_by == "digestions":
                if enzymes_label not in all_patterns:
                    all_patterns[enzymes_label] = OrderedDict()
                all_patterns[enzymes_label][record.id] = (bands, cuts)
            else:
                if record.id not in all_patterns:
                    all_patterns[record.id] = OrderedDict()
                all_patterns[record.id][enzymes_label] = (bands, cuts)

    Y = len(all_patterns)
    X = len(list(all_patterns.values())[0])
    if axes is None:
        _, axes = plt.subplots(Y, 1, figsize=(0.9 * X, 3 * Y))
    if Y == 1:
        axes = [axes]
    bands_props = {"band_thickness": 2.5}
    if show_band_sizes:
        bands_props.update(dict(label="=size", label_fontdict=dict(size=6)))
    for ax, (cat1, cat2s) in zip(axes, sorted(all_patterns.items())):
        patterns = [
            BandsPattern(
                _bands,
                ladder=ladder,
                label=cat2 if (ax == axes[0]) else None,
                label_fontdict=dict(rotation=70),
                global_bands_props=bands_props,
                band_is_uncut=len(_cuts) == 0,
            )
            for cat2, (_bands, _cuts) in cat2s.items()
        ]
        if plot_ladder:
            patterns = [ladder] + patterns
        pattern_set = BandsPatternsSet(
            patterns=patterns,
            ladder=ladder,
            ladder_ticks=4,
            ticks_fontdict=dict(size=9),
            label=cat1,
        )
        pattern_set.plot(ax)

    return axes
