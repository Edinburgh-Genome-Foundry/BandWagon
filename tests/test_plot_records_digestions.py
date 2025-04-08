"""Basic run-through tests for plotting functions."""

from bandwagon import (
    plot_all_digestion_patterns,
    plot_records_digestions,
    load_record,
    LADDER_100_to_4k,
)
import os

records_folder = os.path.join("tests", "data", "records")
records = [
    load_record(
        os.path.join(records_folder, filename), id=filename, topology="circular"
    )
    for filename in os.listdir(records_folder)
]
digestions = [("BamHI", "NcoI"), ("BsaI", "XbaI"), ("StyI",)]


def test_plot_records_digestions(tmpdir):
    target = os.path.join(str(tmpdir), "test.pdf")
    plot_records_digestions(
        records=records, digestions=digestions, ladder=LADDER_100_to_4k, target=target
    )


def test_plot_records_digestions_ziplist(tmpdir):
    target = os.path.join(str(tmpdir), "test.pdf")
    plot_records_digestions(
        records_and_digestions=list(zip(records, digestions)),
        ladder=LADDER_100_to_4k,
        target=target,
    )


def test_plot_all_digestion_patterns():
    plot_all_digestion_patterns(
        records=records, digestions=digestions, ladder=LADDER_100_to_4k
    )
