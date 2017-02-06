import os
import itertools
import pytest
import matplotlib
matplotlib.use("Agg")

from bandwagon import (compute_digestion_bands, BandsPattern,
                       BandsPatternsSet, LADDER_100_to_4k)

with open(os.path.join("tests", "example_sequence.txt"), "r") as f:
    example_sequence = f.read()
ladder = LADDER_100_to_4k.modified(label="Ladder", background_color="#ffebe6")

@pytest.mark.parametrize("digestion,linear,expected", [
    (["EcoRI"], True, [400, 1017, 3583]),
    (["EcoRI"], False, [1017, 3983]),
    (["EcoRI", "BamHI"], True, [400, 417, 600, 3583])
])
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
    ]
    patterns_set = BandsPatternsSet(patterns=[ladder] + patterns,
                                    ladder=ladder,
                                    label="Digests", ladder_ticks=3)
    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)


def test_simple_digestions(tmpdir):
    png_path = os.path.join(str(tmpdir), "test.png")
    patterns = [
        BandsPattern(compute_digestion_bands(example_sequence, [enzyme],
                                             linear=True),
                     ladder=LADDER_100_to_4k, label=enzyme)
        for enzyme in ["BamHI", "EcoRI", "EcoRV", "PstI", "SpeI", "XbaI"]
    ]
    patterns_set = BandsPatternsSet(patterns=[LADDER_100_to_4k] + patterns,
                                    ladder=LADDER_100_to_4k,
                                    label="Digestion results", ladder_ticks=3)

    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)


def test_mixed_digestions(tmpdir):
    png_path = os.path.join(str(tmpdir), "test.png")
    enzymes = "EcoRI", "EcoRV", "BamHI", "XhoI"
    mixes = [[e] for e in enzymes] + list(itertools.combinations(enzymes, 2))

    patterns = [
        BandsPattern(compute_digestion_bands(example_sequence, mix,
                                             linear=True),
                     ladder=LADDER_100_to_4k, label=" + ".join(mix),
                     label_fontdict={"rotation": 40, "size": 9})
        for mix in mixes
    ]
    patterns_set = BandsPatternsSet(patterns=[LADDER_100_to_4k] + patterns,
                                    ladder=LADDER_100_to_4k,
                                    label="Digestion results", ladder_ticks=3)
    ax = patterns_set.plot()
    ax.figure.savefig(png_path, bbox_inches="tight", dpi=200)
    assert os.path.exists(png_path)
