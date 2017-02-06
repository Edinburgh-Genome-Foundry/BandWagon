"""Simple gel simulation with gelsimulation.py.

This example shows how to plot the digestion patterns produced by different
restriction enzymes on a same DNA sequence.
"""

from bandwagon import (BandsPattern, BandsPatternsSet, LADDER_100_to_4k,
                       compute_digestion_bands)

with open("example_sequence.txt", "r") as f:
    sequence = f.read()

patterns = [
    BandsPattern(compute_digestion_bands(sequence, [enzyme], linear=True),
                 ladder=LADDER_100_to_4k, label=enzyme)
    for enzyme in ["BamHI", "EcoRI", "EcoRV", "PstI", "SpeI", "XbaI"]
]
patterns_set = BandsPatternsSet(patterns=[LADDER_100_to_4k] + patterns,
                                ladder=LADDER_100_to_4k,
                                label="Digestion results", ladder_ticks=3)

ax = patterns_set.plot()
ax.figure.savefig("simple_digestions.png", bbox_inches="tight", dpi=120)
