import os
import itertools
import pytest
import matplotlib
from bokeh.resources import CDN
from bokeh.embed import file_html

matplotlib.use("Agg")
import matplotlib.image as mpimg

import Bio
from bandwagon import (
    compute_digestion_bands,
    BandsPattern,
    BandsPatternsSet,
    LADDER_100_to_4k,
)
from bandwagon.tools import (
    random_dna_sequence,
    sequence_to_biopython_record,
    load_record,
    record_is_linear,
)

with open(os.path.join("tests", "data", "example_sequence.txt"), "r") as f:
    example_sequence = f.read()
ladder = LADDER_100_to_4k.modified(label="Ladder", background_color="#ffebe6")


def test_random_dna_sequence():
    assert len(random_dna_sequence(5)) == 5
    assert len(set(random_dna_sequence(15)) | set(["A", "T", "C", "G"])) == 4


def test_sequence_to_biopython_record():
    record = sequence_to_biopython_record("ACTAGTCGATGAC")
    assert type(record) == Bio.SeqRecord.SeqRecord
    assert record.seq == "ACTAGTCGATGAC"


def test_load_record():
    with pytest.raises(ValueError):
        assert load_record("test.py")
    record = load_record(
        os.path.join("tests", "data", "records", "asm_00.gb"),
        file_format="gb",
        topology="auto",
        id="auto",
    )
    assert type(record) == Bio.SeqRecord.SeqRecord
    assert record.annotations["topology"] == "linear"
    assert record.id == "asm_00"

    del record.annotations["topology"]
    assert not record_is_linear(record, default=False)


@pytest.mark.parametrize(
    "digestion,linear,expected",
    [
        (["EcoRI"], True, [400, 1017, 3583]),
        (["EcoRI"], False, [1017, 3983]),
        (["EcoRI", "BamHI"], True, [400, 417, 600, 3583]),
    ],
)
def test_compute_digestion_bands(digestion, linear, expected):
    bands = compute_digestion_bands(example_sequence, digestion, linear=linear)
    assert set(expected) == set(bands)


def test_simple_band_pattern(tmpdir):
    png_path = os.path.join(str(tmpdir), "test.png")
    patterns = [
        BandsPattern([100, 500, 3500], ladder, label="C1"),
        BandsPattern([300, 400, 1500], ladder, label="C2"),
        BandsPattern([100, 1200, 1400, 3000], ladder, label="C3"),
        BandsPattern([100, 700], ladder, label="C4"),
        BandsPattern([200], ladder, label="C4", band_is_uncut=True),
    ]
    patterns_set = BandsPatternsSet(
        patterns=[ladder] + patterns, ladder=ladder, label="Digests", ladder_ticks=3
    )
    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)


def test_simple_digestions_matplotlib_plot(tmpdir):
    png_path = os.path.join(str(tmpdir), "test.png")
    patterns = [
        BandsPattern(
            compute_digestion_bands(example_sequence, [enzyme]),
            ladder=LADDER_100_to_4k,
            label=enzyme,
            global_bands_props={"label": "=size"},
        )
        for enzyme in ["BamHI", "EcoRI", "EcoRV", "PstI", "SpeI", "XbaI"]
    ]
    patterns_set = BandsPatternsSet(
        patterns=[LADDER_100_to_4k] + patterns,
        ladder=LADDER_100_to_4k,
        label="Digestion results",
        ladder_ticks=3,
    )

    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)


def test_simple_digestions_bokeh_plot(tmpdir):
    patterns = [
        BandsPattern(
            compute_digestion_bands(example_sequence, [enzyme]),
            ladder=LADDER_100_to_4k,
            label=enzyme,
            global_bands_props={"label": "=size"},
        )
        for enzyme in ["BamHI", "EcoRI", "EcoRV", "PstI", "SpeI", "XbaI"]
    ]
    patterns_set = BandsPatternsSet(
        patterns=[LADDER_100_to_4k] + patterns,
        ladder=LADDER_100_to_4k,
        label="Digestion results",
        ladder_ticks=3,
    )

    plot = patterns_set.plot_with_bokeh()
    target_file = os.path.join(str(tmpdir), "plot_with_bokeh.html")
    with open(target_file, "w+") as f:
        f.write(file_html(plot, CDN, "Example Sequence"))
    with open(target_file, "r") as f:
        assert len(f.read()) > 9400


def test_mixed_digestions(tmpdir):
    png_path = os.path.join(str(tmpdir), "test.png")
    enzymes = "EcoRI", "EcoRV", "BamHI", "XhoI"
    mixes = [[e] for e in enzymes] + list(itertools.combinations(enzymes, 2))

    patterns = [
        BandsPattern(
            compute_digestion_bands(example_sequence, mix, linear=True),
            ladder=LADDER_100_to_4k,
            label=" + ".join(mix),
            label_fontdict={"rotation": 40, "size": 9},
        )
        for mix in mixes
    ]
    patterns_set = BandsPatternsSet(
        patterns=[LADDER_100_to_4k] + patterns,
        ladder=LADDER_100_to_4k,
        label="Digestion results",
        ladder_ticks=3,
    )
    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)


def test_customized_band_pattern(tmpdir):
    band_image_path = os.path.join("tests", "data", "band_image.png")
    png_path = os.path.join(str(tmpdir), "test.png")
    gel_image = mpimg.imread(band_image_path)
    patterns = [
        BandsPattern(
            [100, 1500, 2000, 1000, 3500],
            ladder,
            label="C1",
            corner_note="corner note",
            corner_note_fontdict=None,
            background_color="#eeeeff",
            gel_image=gel_image,
        )
    ]
    patterns_set = BandsPatternsSet(
        patterns=[ladder] + patterns, ladder=ladder, label="Digests", ladder_ticks=3
    )
    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)
