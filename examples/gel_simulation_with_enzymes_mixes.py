"""Gel simulation with emzymatic mixes and custom ladder.

This example extends `simple_gel_simulation.py` with mixes of enzymes and
a user-defined ladder obtained by measuring the y-coordinate of the ladder's
bands in a gel picture.
"""

from gelsimulator import GelSimulator, GelLadder
import matplotlib.pyplot as plt
import itertools

with open("example_sequence.txt", "r") as f:
    sequence = f.read()

enzymes_names = ["EcoRI", "EcoRV", "BamHI"]
enzymes_mixes = ([[e] for e in enzymes_names] +
                 list(itertools.combinations(enzymes_names, 2)))

custom_ladder = GelLadder(bands={
    # band_size : meaured y-coordinate
    100: 200,
    300: 170,
    500: 150,
    1650: 100,
    4000: 65
})
gel_simulator = GelSimulator(custom_ladder)


fig, ax = plt.subplots(1, figsize=(1 + 1.2 * len(enzymes_mixes), 5))
gel_simulator.format_ax(ax)
gel_simulator.plot_ladder(ax, x_coord=1)
for i, enzymes in enumerate(enzymes_mixes):
    label = " +\n".join(enzymes)
    gel_simulator.plot_digestion_result(ax, sequence, enzymes,
                                        x_coord=i + 2, label=label)
fig.tight_layout()
fig.savefig("multi_enzyme_mixes.png", bbox_inches="tight")
